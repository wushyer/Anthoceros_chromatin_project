#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=m
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --time=3-00:00:00
#SBATCH --output=my.stdout.hicexp1

module load anaconda3/2019.03
conda activate /software/extra/hicexplorer
module load bwa/0.7.17-foss-2018b
module load samtools/1.9-foss-2018b

bwa index Aa.us.chr.fa
bwa mem -A 1 -B 4 -E 50 -L 0 -t 8 Aa.us.chr.fa Anthoceros_agrestis_OXF_contamfiltered_R1.fastq.gz  | samtools view -Shb - > R1.bam
bwa mem -A 1 -B 4 -E 50 -L 0 -t 8 Aa.us.chr.fa Anthoceros_agrestis_OXF_contamfiltered_R2.fastq.gz | samtools view -Shb - > R2.bam


hicBuildMatrix --samFiles R1.bam R2.bam --binSize 10000 --restrictionSequence GATC --danglingSequence GATC --outBam hic_R1.bam --outFileName hic_R1_10kb.h5 --QCfolder hic_R1_10kb_QC --inputBufferSize 5000000


hicMergeMatrixBins --matrix hic_R1_10kb.h5 --numBins 2 --outFileName hic_R1.20kb.h5
hicMergeMatrixBins --matrix hic_R1_10kb.h5 --numBins 4 --outFileName hic_R1.40kb.h5
hicMergeMatrixBins --matrix hic_R1_10kb.h5 --numBins 10 --outFileName hic_R1.100kb.h5
hicMergeMatrixBins --matrix hic_R1_10kb.h5 --numBins 20 --outFileName hic_R1.200kb.h5

hicCorrectMatrix diagnostic_plot --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 -m hic_R1.20kb.h5 -o hic_R1.20kb.png
hicCorrectMatrix diagnostic_plot --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 -m hic_R1.40kb.h5 -o hic_R1.40kb.png
hicCorrectMatrix diagnostic_plot --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 -m hic_R1.100kb.h5 -o hic_R1.100kb.png
hicCorrectMatrix diagnostic_plot --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 -m hic_R1.200kb.h5 -o hic_R1.200kb.png


hicCorrectMatrix correct --chromosomes chr1 chr2 chr3 chr4 chr5 chr6 --filterThreshold -2 2 --perchr -m hic_R1.40kb.h5 --sequencedCountCutoff 0.2 --iterNum 1500 -o hicexplorer.40kb.Corrected.h5
hicFindTADs -m hicexplorer.40kb.Corrected.h5 --correctForMultipleTesting fdr --delta 0.01 --outPrefix hicexplorer.40kb.Corrected-default --numberOfProcessors 8
hicFindTADs -m hicexplorer.40kb.Corrected.h5 --minDepth 200000 --maxDepth 3000000 --step 40000 --minBoundaryDistance 400000 --correctForMultipleTesting fdr --delta 0.01 --outPrefix hicexplorer.40kb.Corrected --numberOfProcessors 8


####note check the last gene region
hicTransform -m ../hicexplorer.40kb.Corrected.h5 --outFileName obs-exp-lib.h5 --method obs_exp_lieberman
hicPCA -m ../hicexplorer.40kb.Corrected.h5 -o pca3.bw pca4.bw --format bigwig --extraTrack gene.bed
hicPlotMatrix -m obs-exp-lib.h5 --outFileName nogene-pca1.pdf --perChr --bigwig pca3.bw --vMin 1 --vMax 10 --dpi 300

