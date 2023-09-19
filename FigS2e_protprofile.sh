#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=logs/%x_%A_%a.txt
# === end SBATCH directives ===
cd /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig
mkdir -p /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/profile_plot
module load deeptools/3.3.1-foss-2018b-python-3.6.6
computeMatrix reference-point \
  --referencePoint TSS \
  -p 8 \
  -b 2000 \
  -a 2000 \
  -bs 10 \
  -R /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS/AagrOXFv2.Gene.cluster.k5.mod.all.bed \
  -S \
	      H3K9me1.log2r_H3_uniq_merged.bw \
        H3K27me1.log2r_H3_uniq_merged.bw \
	      H3K4me3.log2r_H3_uniq_merged.bw \
	      H3K36me3.log2r_H3_uniq_merged.bw \
        H3K27me3.log2r_H3_uniq_merged.bw \
        --skipZeros \
	-o profile_plot/reference_point_TSS_2000.genes.gz 

  plotProfile -m profile_plot/reference_point_TSS_2000.genes.gz \
         -out profile_plot/reference_point_TSS_2000.genes.profile.pdf \
         --outFileSortedRegions profile_plot/reference_point_TSS_2000.genes.profile.bed  \
         --perGroup \
         --colors '#882255' '#CC6677' '#117733' '#44AA99' '#332288' \
         --samplesLabel 'H3K9me1' 'H3K27me1' 'H3K4me3' 'H3K36me3' 'H3K27me3' 
         