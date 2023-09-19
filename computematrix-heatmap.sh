#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=%A_%a.txt
# === end SBATCH directives ===


# === begin ENVIRONMENT SETUP ===

#4. set directory containing sample folders
work_folder=/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig
# Load the required modules
module load deeptools/3.3.1-foss-2018b-python-3.6.6
# module load deeptools/3.1.2-foss-2018b-python-2.7.15

# ... and then change to  working directory:
cd $work_folder
make directory for outputs
mkdir -p $work_folder/reference_point_TSS

computeMatrix reference-point \
  --referencePoint TSS \
  -p 8 \
  -b 2000 \
  -a 2000 \
  -bs 10 \
  -R /groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2/AagrOXF.gene.v2.bed \
  -S \
	H3K9me1.log2r_H3_uniq_merged.bw \
    H3K27me1.log2r_H3_uniq_merged.bw \
	H3K4me3.log2r_H3_uniq_merged.bw \
	H3K36me3.log2r_H3_uniq_merged.bw \
	H3K27me3.log2r_H3_uniq_merged.bw \
  --skipZeros \
	-o reference_point_TSS/reference_point_TSS_2000.genes.gz \
 
for i in 3 4 5 6 7 8
 do  
 plotHeatmap -m reference_point_TSS/reference_point_TSS_2000.genes.gz \
         -out reference_point_TSS/reference_point_TSS_2000.genes.k$i.pdf \
         --colorMap coolwarm \
         --zMin -2 --zMax 2 \
         --kmeans $i \
         --samplesLabel H3K9me1 H3K27me1 H3K4me3 H3K36me3  H3K27me3 \
         --outFileSortedRegions reference_point_TSS/reference_point_TSS_2000.genes.k$i.bed
done

computeMatrix scale-regions \
  -p 8 \
  --regionBodyLength 1000 \
  -b 1000 \
  -a 1000 \
  -bs 10 \
  -R /groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2/Aa.us.TE.bed \
  -S \
	H3K9me1.log2r_H3_uniq_merged.bw \
  H3K27me1.log2r_H3_uniq_merged.bw \
	H3K4me3.log2r_H3_uniq_merged.bw \
	H3K36me3.log2r_H3_uniq_merged.bw \
	H3K27me3.log2r_H3_uniq_merged.bw \
  --skipZeros \
	-o scale_regions/scale_regions1000.TEs.gz \
 
 for i in 3 4 5 6 7 8
 do  
 plotHeatmap -m scale_regions/scale_regions1000.TEs.gz \
         -out scale_regions/scale_regions1000.TEs.heatmap.k$i.pdf \
         --colorMap coolwarm \
         --zMin -2 --zMax 2 \
         --kmeans $i \
         --samplesLabel H3K9me1 H3K27me1 H3K4me3 H3K36me3  H3K27me3 \
         --outFileSortedRegions scale_regions/scale_regions1000.TEs.heatmap.k$i.bed
done
