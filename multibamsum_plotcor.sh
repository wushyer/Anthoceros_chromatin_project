#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --output=%A_%a.txt
#SBATCH --error=%A_%a.err
# === end SBATCH directives ===


cd /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq
module load deeptools/3.3.1-foss-2018b-python-3.6.6
multiBamSummary bins --bamfiles */*.sorted_uniq.bam -out readCounts.np --smartLabels

plotCorrelation \
    -in readCounts.np \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_SpearmanCorr_readCounts.png   \
    --outFileCorMatrix SpearmanCorr_readCounts.tab

plotCorrelation \
    -in readCounts.np \
    --corMethod pearson --skipZeros \
    --plotTitle "Pearson Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_PeasonCorr_readCounts.png   \
    --outFileCorMatrix PeasonCorr_readCounts.tab
