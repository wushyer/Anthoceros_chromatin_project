#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=logs/%x_%A_%a.txt
#SBATCH --error=logs/%x_%A_%a.err
# === end SBATCH directives ===

#Script to trim adapters then align CUT&RUN paired end data with Bowtie2
#Script adapted from Michael Borg

# === begin ENVIRONMENT SETUP ===

#1. list of sample names on separate lines (eg copied from excel)
#sample_list=/groups/berger/lab/tetsuya.hisanaga/input_files/EMseq_AaOXF_ga.txt
#5. Genome mapped to
output=AagrOXFv2 #size=208291648
chromsize=127252112
# output=vector_spikein
#4. set directory containing sample folders
work_folder=/scratch-cbe/users/tetsuya.hisanaga/EMseq/${output}
#5. bwt2 index for genome of interest
index1=/groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/${output}
#6. Data type "SE" or "PE"
datatype="PE"

#F=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
#4. sample name
#NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $2}'`

# Load the required modules
module load bismark/0.22.2-foss-2018b
module load samtools/1.9-foss-2018b
module load bedtools/2.27.1-foss-2018b
module load cutadapt/1.18-foss-2018b-python-3.6.6
module load bowtie2/2.3.4.2-foss-2018b
module load deeptools/2.5.4-foss-2018b-python-2.7.15
module load picard/2.18.27-java-1.8
module load fastqc/0.11.8-java-1.8
module load trim_galore/0.6.2-foss-2018b-python-3.6.6

# === end ENVIRONMENT SETUP ===
mkdir -p $work_folder/merged
cd $work_folder

export TMPDIR=$work_folder/tmp
# make sure the directory exists
mkdir -p $TMPDIR

##########################################################################################

samtools merge merged/AaOXF_ga_merged.bam AaOXF_ga_rep3/HGYF2DRXY_2#161179_TACTCCAGGGTATAGG.bam AaOXF_ga_rep8/HGYF2DRXY_2#161180_AGTGACCTTCTAGGAG.bam
samtools sort -n -o merged/AaOXF_ga_merged.qsort -O bam merged/AaOXF_ga_merged.bam

NAME=AaOXF_ga_merged
cd $work_folder/merged

##Convert bam to fastq
 if [ $datatype == "PE" ]; then
         bedtools bamtofastq -i ${NAME}.qsort -fq ${NAME}.end1.fq -fq2 ${NAME}.end2.fq
 fi
 if [ $datatype == "SE" ]; then
        bedtools bamtofastq -i ${NAME}.qsort -fq ${NAME}.fq
 fi
## Remove adaptors
 trim_galore --illumina --paired --clip_R2 2 ${NAME}.end1.fq ${NAME}.end2.fq -o $work_folder/merged --gzip

 ## Prepare genome for bismark
 bismark_genome_preparation ${index1}

 ## Run bismark
 bismark --genome ${index1} --bowtie2 --seedmms 1 -1 ${NAME}.end1_val_1.fq.gz -2 ${NAME}.end2_val_2.fq.gz --fastq --multicore 4
 ##Deduplicate bam
 deduplicate_bismark --bam -p ${NAME}.end1_val_1_bismark_bt2_pe.bam
#Extract methylation information
 bismark_methylation_extractor --cytosine_report --CX --bedGraph --multicore 4 --gzip  --comprehensive --no_header --genome_folder ${index1} -p ${NAME}.end1_val_1_bismark_bt2_pe.deduplicated.bam 
#Create report
 bismark2report
 #bismark2summary
