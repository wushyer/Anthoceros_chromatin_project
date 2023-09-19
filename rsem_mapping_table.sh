#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/script/20230124AaRNAseq_mapping_table_1.txt
# === end SBATCH directives ===

#Adapted from Nodine RSEM pipeline
#load samtools module and redirect out put to null to give a clean ouput
ml samtools/1.9-foss-2018b &> /dev/null

#1. set directory containing sample folders
input_dir=/scratch-cbe/users/tetsuya.hisanaga/20230124AaRNAseq1/ #folder containing subfolders of samples
x=1 #counter later used for pbs array index number

#TODO implement some system that autodetects the adaptor_type
adaptor_type="unknown" #default adaptor_type
  #loop over subfolders and identify the seq file types in the subfolders
for d in $(readlink -m $input_dir/*/) ;
do
  #initilize some variables
  file_type=""
  seq_type=""
  check_for_bam="no"

  fq=($(ls $d | grep -e ".fq.gz\|.fastq.gz\|.fq\|.fastq")) #grep all fq files
  #some sanity checking...
  if [[ "${#fq[@]}" -gt "2" ]]; then
    echo "Warning! More than 2 fq files found! Skipping... Check path $d" >&2
  #when 2 fq files are found mode is PE and fq
  elif [[ "${#fq[@]}" -eq "2" ]]; then
    seq_type="PE"
    file_type="fq"
  #when 1 fq file is found mode is SE and fq
  elif [[ "${#fq[@]}" -eq "1" ]]; then
    seq_type="SE"
    file_type="fq"
  #when no fq files are found check for bam files
  else
    check_for_bam="yes"
  fi

  if [[ "$check_for_bam" == "yes" ]]; then
    bam=($(ls $d | grep -e ".bam")) #grep all bam files
    #some sanity checking...
    if [[ "${#bam[@]}" -gt "1" ]]; then
      echo "Warning! More than 1 bam file found! Skipping... Check path $d" >&2
    #when a bam file is found set the mode to bam end test the seq_type
    elif [[ "${#bam[@]}" -eq "1" ]]; then
      file_type="bam"
      #check with samtools flag how many paired reads are present
      pe_reads=$(samtools view -c -f 1 $d/$bam)
      #when paired reads are present set mode to PE, otherwise set it to SE
      if [[ $pe_reads -gt "0" ]]; then
        seq_type="PE"
      else
        seq_type="SE"
      fi
    else
      #print some error when neither fq nor bam files are found
      echo "Warning! Neither fq nor bam files found! Skipping... Check path $d" >&2
    fi
  fi

  printf "%d %s %s %s %s\n" $x $d $file_type $seq_type $adaptor_type 
  ((x++))
done
