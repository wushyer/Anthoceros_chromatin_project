#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=8
#SBATCH --output=%A_%a.txt
#SBATCH --error=%A_%a.err
# === end SBATCH directives ===

# === begin ENVIRONMENT SETUP ===

#1. list of sample names on separate lines (eg copied from excel)
sample_list=/groups/berger/lab/tetsuya.hisanaga/input_files/20230123ChIP.txt

#4. set directory containing sample folders
work_folder=/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq
#6. line in sample file input sample name is on

#4. sample name
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`


# Load the required modules
module load samtools/1.9-foss-2018b

# ... and then change to  working directory:
cd $work_folder
# === end ENVIRONMENT SETUP ===


#set variables
mkdir merged_bam
#merge replicates
samtools merge merged_bam/${NAME}.1.sorted_uniq_merged.bam ${NAME}_rep3/${NAME}_rep3.sorted_uniq.bam ${NAME}_rep4/${NAME}_rep4.sorted_uniq.bam
samtools index merged_bam/${NAME}.1.sorted_uniq_merged.bam
