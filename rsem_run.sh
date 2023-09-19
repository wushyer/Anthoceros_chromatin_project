#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1-12
#SBATCH --output=logs/%x_%A_%a.txt
# === end SBATCH directives ===

#TODO: update README & make the make_pbs_mapping file.
#Adapted from Nodine RSEM pipeline

# === begin ENVIRONMENT SETUP ===
####set to 0 (false) or 1 (true) to let the repsective code block run
#1. run rsem
#TODO: add conditional adaptor trimming
#1. trim adaptors
trim_adaptors=1
#2. run rsem
run_rsem=1
#3. make plots or not
make_plots=1
#4. delete unecessary files from temp_dir
clean=1
##### specify RSEM parameters
aligner=star
threads=4 #set this to the number of available cores
##### specify folders and variables #####
#set script dir
pipe_dir=/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/script
#set ouput base dir -> do not need to change
base_dir=/scratch-cbe/users/tetsuya.hisanaga/20230124AaRNAseq1/
#folder for rsem reference
rsem_ref_dir=/groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2/AagrOXF_rsem-index
#add folder basename as prefix (follows convention from rsem_make_reference) -> do not change
rsem_ref=$rsem_ref_dir/$(basename "$rsem_ref_dir")
#location of the mapping file for the array job -> make sure it is the same as output from rsem_mapping_table.sh
pbs_mapping_file=/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/script/20230124AaRNAseq_mapping_table_1.txt
#super folder of the temp dir, script will create subfolders with $sample_name -> do not change
temp_dir=$base_dir/tmp


if [ $make_plots -eq 1 ]; then
  module load r/3.5.1-foss-2018b
fi
##### Obtain Parameters from mapping file using $PBS_ARRAY_INDEX as line number
input_mapper=`sed -n "${SLURM_ARRAY_TASK_ID} p" $pbs_mapping_file` #read mapping file
names_mapped=($input_mapper)

sample_dir=${names_mapped[1]} # get the sample dir
file_type=${names_mapped[2]} # get the file type
seq_type=${names_mapped[3]} # get the seq type
adaptor_type=${names_mapped[4]} # get the adaptor type

sample_name=`basename $sample_dir` #get the base name of the dir as sample name

#print some output for logging
echo '#########################################################################'
echo 'Starting RSEM RNA-seq pipeline for: ' $sample_name
echo 'Sample directory: ' $sample_dir
echo 'Rsem reference: ' $rsem_ref
echo 'Aligner to be used: ' $aligner
echo 'Mapping file: ' $pbs_mapping_file
echo 'Specified file type: ' $file_type
echo 'Specified sequencing mode: ' $seq_type
echo 'Specified adaptor type: ' $adaptor_type
echo '#########################################################################'

#some error handling function
function error_exit
{
  echo "$1" 1>&2
  exit 1
}

#make output folder
cd $sample_dir

#folders for temp files
temp_dir_s=$temp_dir/$sample_name
mkdir -p $temp_dir_s

samples_trimmed=$sample_dir/trimmed
mkdir -p $samples_trimmed

if [ $run_rsem -eq 1 ]; then
  #1. check file typ and convert to fastq
  case $file_type in
    "bam")
      echo "File type bam. Converting to fastq..."
      f=($(ls $sample_dir | grep -e ".bam")) # get all bam files in folder
      #throw error if more or less than 1 file is present
      if [[ "${#f[@]}" -ne "1" ]]; then
        error_exit "Error: wrong number of bam files in folder. Files present: ${#f[@]}"
      fi
      #load modules
      module load samtools/1.9-foss-2018b
      module load bedtools/2.27.1-foss-2018b

      samtools sort -n -m 4G -@ $threads -o $sample_dir/${f%.*}.sorted.bam \
        $sample_dir/$f
      bedtools_params="bedtools bamtofastq -i $sample_dir/${f%.*}.sorted.bam "
      case $seq_type in
        "PE") #modify bedtools params for PE conversion
          bedtools_params=$bedtools_params" -fq $sample_dir/${f%.*}.1.fq"\
            " -fq2 $sample_dir/${f%.*}.2.fq"
          ;;
        "SE") #modify bedtools params for SE conversion
          bedtools_params=$bedtools_params" -fq $sample_dir/${f%.*}.fq"
          ;;
        *) #exit when unexpected input is encountered
          error_exit "Error: wrong paramter for seq type selected! Select PE or SE."
          ;;
      esac
      #print the command to be exectuted
      echo "Command exectuted for converting bam to fastq:\n $bedtools_params"
      eval $bedtools_params #run the command
      echo "Converting to fastq... Done"
      ;;
    "fq")
      echo "File type fastq. No conversion necessary..."
      # do nothing...
      ;;
    *) #exit when unexpected input is encountered
      error_exit "Error: wrong paramter for file type selected! Select bam or fq."
      ;;
  esac

  #2. do adaptor trimming according to seq_type and adaptor_type
  #get the zipped files
  f=($(ls $sample_dir | grep -e ".fq.gz\|.fastq.gz"))
  #check if more than 0 zipped files are present, if so unzip
  if [[ "${#f[@]}" -gt "0" ]]; then
    gunzip "${f[@]}"
  fi

  f=($(ls $sample_dir| grep -e ".fq\|.fastq"))
  trim_params="trim_galore --dont_gzip --stringency 4 -o $samples_trimmed"
  case $adaptor_type in
    "nextera")
      trimming=$trimming" --nextera"
      ;;
    "illumina")
      trimming=$trimming" --illumina"
      ;;
    "unknown") #run trim galore with autodetect
      #do nothing == autodetect
      ;;
    "none")
      #FIXME this will cause the pipeline to crash atm.
      trim_params="No trimming selected..." #Don't trimm
      ;;
    ^[NCAGTncagt]+$) #check if alphabet corresponds to the genetic alphabet
      if [[ $seq_type == "SE" ]]; then
        trimming=$trimming" -a $adaptor_type"
      else
        error_exit "Error: Wrong paramter for adaptor or seq type selected! \
          See documentation for valid types"
      fi
      ;;
    ^[NCAGTncagt\/]+$) #check if alphabet corresponds to the genetic alphabet
      if [[ $seq_type == "PE" ]]; then
        seqs=(${adaptor_type//\// })
        trimming=$trimming" -a ${seqs[0]} -a2 ${seqs[1]}"
      else
        error_exit "Error: Wrong paramter for adaptor or seq type selected! \
          See documentation for valid types"
      fi
      ;;
    *) #exit when unexpected input is encountered
      error_exit "Error: Wrong paramter for adaptor type selected! \
        See documentation for valid types"
      ;;
  esac

  fq1=""
  fq2=""
  fq=""
  case $seq_type in
    "PE")
      trim_params=$trim_params" --paired $sample_dir/${f[0]} $sample_dir/${f[1]}"
      #fq1=$sample_dir/$samples_trimmed/${f[0]%.*}1_val_1.fq
      #fq2=$sample_dir/$samples_trimmed/${f[1]%.*}2_val_2.fq
      fq1=$samples_trimmed/${f[0]%.*}_val_1.fq
      fq2=$samples_trimmed/${f[1]%.*}_val_2.fq
      ;;
    "SE")
      trim_params=$trim_params" $sample_dir/$f"
      #TODO: check if the file name is correct like this.
      #fq=$sample_dir/$samples_trimmed/${f%.*}_trimmed.fq
      fq=$samples_trimmed/${f%.*}_trimmed.fq
      ;;
    *) #exit when unexpected input is encountered
      error_exit "Error: Wrong paramter for seq type selected! \
        See documentation for valid types"
      ;;
  esac
  #load module
  module load trim_galore/0.6.2-foss-2018b-python-3.6.6
  #print the command to be exectuted
  echo "Command exectuted for adaptor trimming:" \n "$trim_params"
  if [[ $adaptor_type != "none" ]]; then
     eval "$trim_params" #run the command
  fi

  rsem_opts=""
  case $seq_type in
    "PE")
      rsem_opts=$rsem_opts"--paired-end $fq1 $fq2"
      ;;
    "SE")
      rsem_opts=$rsem_opts"$fq"
      ;;
    *) #exit when unexpected input is encountered
      error_exit "Error: Wrong paramter for seq type selected!" \
        "See documentation for valid types"
      ;;
  esac

  # run rsem to calculate the expression levels
  # --estimate-rspd: estimate read start position to check if the data has bias
  # --output-genome-bam: output bam file as genomic, not transcript coordinates
  # --seed 12345 set seed for reproducibility of rng
  # --calc-ci calcutates 95% confidence interval of the expression values
  # --ci-memory 30000 set memory
  rsem_params="--$aligner \
  --num-threads $threads \
  --temporary-folder $temp_dir_s \
  --append-names \
  --estimate-rspd \
  --output-genome-bam \
  --seed 12345 \
  --calc-ci \
  --ci-memory 40000 \
  --sort-bam-by-coordinate" \
  # --alignEndsType EndToEnd \
  # --outSAMattributes NH HI NM MD \
  # --sort-bam-by-read-name"
  ##remove 4th last for SNPsplit; add last 3 ##Doesn't work for SNPsplit -> cannot output MD flag

  rsem_params+=" "$rsem_opts
  rsem_params+=" "$rsem_ref
  rsem_params+=" "$sample_name

  #sem_cmds=$rsem_params $rsem_opts $rsem_ref $sample_name

  
  mkdir -p $sample_dir/rsem/
  cd $sample_dir/rsem/
  # conditional loading of modules based on aligner to be used by RSEM
  if [ "$aligner" == "bowtie" ]; then
    module load bowtie/1.2.2-foss-2018b
  fi
  if [ "$aligner" == "bowtie2" ]; then
    module load bowtie2/2.3.4.2-foss-2018b
  fi
  if [ "$aligner" == "star" ]; then
    module load star/2.7.1a-foss-2018b
  fi
  module load rsem/1.3.2-foss-2018b
  #rsem command that should be run
  echo "rsem-calculate-expression $rsem_params >& $sample_name.log"
  eval "rsem-calculate-expression $rsem_params >& $sample_name.log"
fi

#run the rsem plot function
if [ $make_plots -eq 1 ]; then
  rsem-plot-model $sample_dir/rsem/$sample_name $sample_dir/rsem/$sample_name.pdf
fi

#delete the temp files
if [ $clean -eq 1 ]; then
  gzip $sample_dir/*.fq $sample_dir/*.fastq
  rm $sample_dir/${f%.*}.sorted.bam
  rm $sample_dir/rsem/*.transcript.bam
  rm -rf $temp_dir_s
  rm -rf $samples_trimmed
fi

echo 'Finished RSEM RNA-seq pipeline for: '$sample_name
