#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=04:00:00
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
#2. specify the line used as control
input_line=1
#4. set directory containing sample folders
work_folder=/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam

# Load the required modules



# ... and then change to  working directory:
cd $work_folder
mkdir macs2
mkdir bigwig



#4. sample name
NAME=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`
INPUT=`sed -n "${input_line} p" $sample_list | awk '{print $1}'`
# call MACS2 broad peaks on the final merged data

module load macs2/2.2.5-foss-2018b-python-3.6.6
# set TMPDIR for this shell session (this will set it to "/lustre/scratch/users/michael.borg/macs_tmp")
export TMPDIR=$work_folder/macs_tmp
# make sure the directory exists
mkdir -p $TMPDIR
# run MACS as usual
macs2 callpeak -t ${NAME}.1.sorted_uniq_merged.bam -c ${INPUT}.1.sorted_uniq_merged.bam\
 --nomodel --nolambda -q 0.01 -f BAM --broad --broad-cutoff 0.1 -g 127252112 -n ${NAME} --outdir macs2/

module load deeptools/3.3.1-foss-2018b-python-3.6.6
#set tmp file for deeptools

export TMPDIR=$work_folder/deeptool_tmp
# make sure the directory exists
mkdir -p $TMPDIR
bamCompare -bs=10 -b1 ${NAME}.1.sorted_uniq_merged.bam -b2 ${INPUT}.1.sorted_uniq_merged.bam --scaleFactorsMethod readCount -o bigwig/${NAME}.log2r_H3_uniq_merged.bw
#bamCompare -bs=5 -b1 H3K27me3.merged.sorted.uniq.bam -b2 H3.merged.sorted.uniq.bam --ratio subtract --scaleFactorsMethod readCount --normalizeTo1x 16728945 -o H3K27me3.subtr.bw



