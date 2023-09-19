#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1-8
#SBATCH --output=/users/tetsuya.hisanaga/logs/DEEP-align-Aa-chip_%A_%a.txt
#SBATCH --error=/users/tetsuya.hisanaga/logs/DEEP-align-Aa-chip_%A_%a.err
# === end SBATCH directives ===

#Script to trim adapters then align ChIP-seq paired end data with Bowtie2
#Script adapted from Michael Borg

# === begin ENVIRONMENT SETUP ===

#1. list of sample names on separate lines (eg copied from excel)
sample_list=/groups/berger/lab/tetsuya.hisanaga/input_files/20220425AaChIPseq_1.txt
#2. set directory containing sample folders
work_folder=/scratch-cbe/users/tetsuya.hisanaga/AaChIPseq
#3. bwt2 index
index1=/groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2/Anthoceros_agrestis_Oxford_genome
#4. Data type "SE" or "PE"
datatype="SE"
#5. taking file names from the sample list
F=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
#6. sample name
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $2}'`


# Load the required modules
module load samtools/1.9-foss-2018b
module load bedtools/2.27.1-foss-2018b
module load cutadapt/1.18-foss-2018b-python-3.6.6
module load bowtie2/2.3.4.2-foss-2018b
module load deeptools/3.3.1-foss-2018b-python-3.6.6
module load picard/2.18.27-java-1.8
module load r/3.5.1-foss-2018b
# ... and then change to  working directory:
cd $work_folder
# === end ENVIRONMENT SETUP ===

mkdir ${NAME}

# BAM file must be sorted/grouped by read name to keep records in the two output FASTQ files in the same order.
# To sort the BAM file by query name use ---> samtools sort -n -o aln.qsort -O bam aln.bam

if [[ ${F##*.} == "bam" ]]; then
	ln -s /groups/berger/lab/Raw/demultiplexed/${F} $work_folder/${NAME}/
	samtools sort -n -o ${NAME}/${NAME}.qsort -O bam ${NAME}/${F}
	if [ $datatype == "PE" ]; then
		bedtools bamtofastq -i ${NAME}/${NAME}.qsort -fq ${NAME}/${NAME}.end1.fq -fq2 ${NAME}/${NAME}.end2.fq
	elif [ $datatype == "SE" ]; then
		bedtools bamtofastq -i ${NAME}/${NAME}.qsort -fq ${NAME}/${NAME}.fastq
	else
		echo "Datatype??"
	fi
	echo "BAM converted to SAM for sample ... "
	echo ${NAME}
	rm ${NAME}/${NAME}.qsort
else
	if [ $datatype == "PE" ]; then
		gunzip ${NAME}/${F}_1.fastq.gz
		gunzip ${NAME}/${F}_2.fastq.gz
		mv ${NAME}/${F}_1.fastq ${NAME}/${NAME}.end1.fq
		mv ${NAME}/${F}_2.fastq ${NAME}/${NAME}.end2.fq
	elif [ $datatype == "SE" ]; then
		gunzip ${NAME}/${F}.fastq.gz
		mv ${NAME}/${F}.fastq ${NAME}/${NAME}.fastq
	else
		echo "Filetype??"
	fi
	echo "I tried my best with this SRA file"
fi

a1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
a2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

# cutadapt -a AACCGGTT -o output.fastq input.fastq --- removing sequences less than 1bp (self-ligated adapters)
echo ".............................................................. "
echo "---> Cutting adapters for sample ... " ${NAME}
if [ $datatype == "PE" ]; then
        cutadapt -a $a1 -A $a2 --minimum-length 1 -o ${NAME}/${NAME}.trim1.fq -p ${NAME}/${NAME}.trim2.fq ${NAME}/${NAME}.end1.fq ${NAME}/${NAME}.end2.fq > ${NAME}/${NAME}_cutadapt.txt
fi
if [ $datatype == "SE" ]; then
        cutadapt -a $a1 --minimum-length 1 -o ${NAME}/${NAME}.trim.fq  ${NAME}/${NAME}.fastq > ${NAME}/${NAME}_cutadapt.txt
        # cutadapt -a $a1 --minimum-length 1 -o ${NAME}/${NAME}.trim.fq  ${NAME}/${NAME}.fastq.gz > ${NAME}/${NAME}_cutadapt.txt
fi

# align with bowtie2 fragments between 100bp and 2000bp with 0 mismatches (-N)
echo ".............................................................. "
echo "---> Aligning sample ... " ${NAME}
if [ $datatype == "PE" ]; then
        bowtie2 --sensitive -p 8 -X 2000 -N 0 --trim5 0 --trim3 0 -x $index1 -1 ${NAME}/${NAME}.trim1.fq -2 ${NAME}/${NAME}.trim2.fq -S ${NAME}/${NAME}.sam  
fi
if [ $datatype == "SE" ]; then
        bowtie2 --sensitive -p 8 -N 0 --trim5 0 --trim3 0 -x $index1 -U ${NAME}/${NAME}.trim.fq -S ${NAME}/${NAME}.sam  
fi

# convert sam to bam
echo ".............................................................. "
echo "---> Converting SAM to BAM for sample ... " ${NAME}
samtools view -bS -o ${NAME}/${NAME}.aligned.bam ${NAME}/${NAME}.sam

#Sorting and indexing aligned BAM file
echo ".............................................................. "
echo "---> Sorting by co-ordinates and indexing sample ... " ${NAME}.aligned.bam 
echo ${NAME}
samtools sort -o ${NAME}/${NAME}.aligned_sorted.bam ${NAME}/${NAME}.aligned.bam
samtools index ${NAME}/${NAME}.aligned_sorted.bam

#use samtools to filter away reads with MAPQ < 10
echo ".............................................................. "
echo "---> Filtering for mapQ > 10m sorting and indexing sample ... " ${NAME}.aligned_sorted.bam 
echo ${NAME}
samtools view -bq 10 ${NAME}/${NAME}.aligned_sorted.bam > ${NAME}/${NAME}.aligned_filtered.bam 
samtools sort -o ${NAME}/${NAME}.sorted.bam ${NAME}/${NAME}.aligned_filtered.bam 
samtools index ${NAME}/${NAME}.sorted.bam

#set tmp file for picard tools
export ptTMPDIR=$work_folder/picard_tmp
# make sure the directory exists
mkdir -p $ptTMPDIR

#removing duplicates and indexing resulting BAM file
echo ".............................................................. "
echo "---> Removing duplicates and indexing to generate ... " ${NAME}.sorted_uniq.bam
java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${NAME}/${NAME}.sorted.bam O=${NAME}/${NAME}.sorted_uniq.bam M=${NAME}/${NAME}.dup_metrics.txt AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
samtools index ${NAME}/${NAME}.sorted_uniq.bam

#set tmp file for deeptools
export TMPDIR=$work_folder/deeptool_tmp
# make sure the directory exists
mkdir -p $TMPDIR

#normalise and convert BAM file to BIGWIG to represent coverage
echo ".............................................................. "
echo "---> Generating BW files for final BAM files with and without duplicates ... "

bamCoverage -b ${NAME}/${NAME}.sorted.bam -o ${NAME}/${NAME}.coverage_dups.bw --normalizeUsing RPGC --effectiveGenomeSize 127252112 --binSize=10 ##AaOXF
bamCoverage -b ${NAME}/${NAME}.sorted_uniq.bam -o ${NAME}/${NAME}.coverage_uniq.bw --normalizeUsing RPGC --effectiveGenomeSize 127252112 --binSize=10 ##AaOXF


# generate alignment summary tables
### for aligned.bam (i,e to get raw mapping information)
echo ".............................................................. "
echo "---> Summarising alignment statistics for sample ... " ${NAME}
samtools stats ${NAME}/${NAME}.aligned_sorted.bam > ${NAME}/${NAME}.aligned_stats.txt
samtools stats ${NAME}/${NAME}.sorted.bam > ${NAME}/${NAME}.sorted_stats.txt
for stat in "raw total sequences:" "reads mapped:" "reads mapped and paired:" "reads properly paired:" "reads duplicated:" "average length:" "maximum length:" "average quality:" "insert size average:" "insert size standard deviation:" "inward oriented pairs:" "outward oriented pairs:" "pairs with other orientation:" "pairs on different chromosomes:" ; do
        grep "^SN" ${NAME}/${NAME}.aligned_stats.txt| grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}/${NAME}.aligned_summary.txt
        grep "^SN" ${NAME}/${NAME}.sorted_stats.txt | grep "$stat" | awk -F $'\t' '{print $2,"\t",$3}' >> ${NAME}/${NAME}.sorted_summary.txt
done
#mtreads=$(samtools idxstats ${NAME}/${NAME}.aligned_sorted.bam | grep "mitochondria" | awk '{print $3}')
#echo -e "reads mapping to mitochondria:\t"$mtreads >> ${NAME}/${NAME}.aligned_summary.txt
#mtreads=$(samtools idxstats ${NAME}/${NAME}.sorted.bam | grep "mitochondria" | awk '{print $3}')
#echo -e "reads mapping to chrM:\t"$mtreads >> ${NAME}/${NAME}.sorted_summary.txt
rm ${NAME}/${NAME}.aligned_stats.txt
rm ${NAME}/${NAME}.sorted_stats.txt


