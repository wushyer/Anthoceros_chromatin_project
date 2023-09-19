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

cd /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/genmap

for i in {1..5}
do 
grep "cluster_$i" AagrOXFv2.Gene.cluster.k5.mod.bed |cut -f1-4 > gcl$i.bed
./bigWigAverageOverBed.html AagrOXFv2_genmap_k50_E0.bw gcl$i.bed genmap_out/gcl$i.txt 
done

for i in {1..8}
do 
grep "cluster_$i" AagrOXFv2.TE.cluster.k8.mod.bed |cut -f1-4 > tcl$i.bed
./bigWigAverageOverBed.html AagrOXFv2_genmap_k50_E0.bw tcl$i.bed genmap_out/tcl$i.txt 
done
