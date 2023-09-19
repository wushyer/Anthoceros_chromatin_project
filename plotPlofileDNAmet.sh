#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=logs/%A_%a_%x.out
#SBATCH --error=logs/%A_%a_%x.out

#created chromatin_peaks.bed using following command
##{ cat H3K9me1_peaks.broadPeak; echo "#H3K9me1 peak"; cat H3K27me1_peaks.broadPeak; echo "#H3K27me1 peak" ; cat H3K4me3_peaks.broadPeak; echo "#H3K4me3 peak" ; cat H3K36me3_peaks.broadPeak; echo "#H3K36me3 peak" ; cat H3K27me3_peaks.broadPeak; echo "#H3K27me3 peak"; } > chromatin_peaks.bed

cd /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/EMseq/merged
module load deeptools/3.3.1-foss-2018b-python-3.6.6

#Plot DNA methylation levels over chromatin modification peaks
computeMatrix reference-point \
--referencePoint center \
--binSize 10 \
-b 1000 -a 1000 \
-S AaOXF_ga_merged_mincov10_CG.bw \
   AaOXF_ga_merged_mincov10_CHG.bw \
   AaOXF_ga_merged_mincov10_CHH.bw  \
-R TE_type.bed \
--skipZeros \
-o merged_all_mincov10_peaks.gz 

plotProfile -m merged_all_mincov10_peaks.gz \
-o merged_all_mincov10_peaks.profile.pdf \
--plotType se --perGroup \
--refPointLabel "Peak center" \
--samplesLabel 'mCG' 'mCHG' 'mCHH' \
--colors 'darkred' 'red' 'tomato' \
--regionsLabel 'H3K9me1 peaks' 'H3K27me1 peaks' 'H3K4me3 peaks' 'H3K36me3 peaks' 'H3K27me3 peaks' \
--numPlotsPerRow 1

#Plot DNA methylation levels over PCG clusters
computeMatrix scale-regions \
--binSize 10 \
--regionBodyLength 2000 -b 1000 -a 1000 \
-S AaOXF_ga_merged_mincov10_CG.bw \
   AaOXF_ga_merged_mincov10_CHG.bw \
   AaOXF_ga_merged_mincov10_CHH.bw  \
-R AagrOXFv2.Gene.cluster.k5.mod.bed \
--skipZeros \
-o merged_all_mincov10_genes.k5.gz 

plotProfile -m merged_all_mincov10_genes.k5.gz \
-o merged_all_mincov10_genes.k5.profile.pdf \
--plotType se --perGroup --startLabel "TSS" --endLabel "TES" \
--samplesLabel 'mCG' 'mCHG' 'mCHH' \
--colors 'darkred' 'red' 'tomato' \
--regionsLabel 'Cluster P1' 'Cluster P2' 'Cluster P3' 'Cluster P4' 'Cluster P5'  \
--numPlotsPerRow 4

#Plot DNA methylation levels over TE clusters
computeMatrix scale-regions \
--binSize 10 \
--regionBodyLength 1000 -b 1000 -a 1000 \
-S AaOXF_ga_merged_mincov10_CG.bw \
   AaOXF_ga_merged_mincov10_CHG.bw \
   AaOXF_ga_merged_mincov10_CHH.bw  \
-R AagrOXFv2.TE.cluster.k8.mod.bed \
-o merged_all_mincov10_TEs.gz 

plotProfile -m merged_all_mincov10_TEs.gz  \
-o merged_all_mincov10_TEs.profile.pdf \
--plotType se --perGroup --startLabel "Start" --endLabel "End" \
--samplesLabel 'mCG' 'mCHG' 'mCHH' \
--colors 'darkred' 'red' 'tomato' \
--regionsLabel 'Cluster T1' 'Cluster T2' 'Cluster T3' 'Cluster T4' 'Cluster T5' 'Cluster T6' 'Cluster T7' 'Cluster T8' \
--numPlotsPerRow 3

#Plot DNA methylation levels over each TE class
computeMatrix scale-regions \
--binSize 10 \
--regionBodyLength 1000 -b 1000 -a 1000 \
-S AaOXF_ga_merged_mincov10_CG.bw \
   AaOXF_ga_merged_mincov10_CHG.bw \
   AaOXF_ga_merged_mincov10_CHH.bw  \
-R TE_typ.3.bed \
--skipZeros \
-o merged_all_mincov10_TEclass.gz 

plotProfile -m merged_all_mincov10_TEclass.gz \
-o merged_all_mincov10_TEclass.profile.pdf \
--plotType se --perGroup \
--startLabel "Start" --endLabel "End" \
--samplesLabel 'mCG' 'mCHG' 'mCHH' \
--colors 'darkred' 'red' 'tomato' \
--numPlotsPerRow 3