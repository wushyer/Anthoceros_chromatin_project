#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --output=logs/%x_%A_%a.txt
#SBATCH --error=logs/%x_%A_%a.txt
# === end SBATCH directives ===
sample_list=/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/macs2_H3/samples
cd /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/macs2_H3
module load  bedtools/2.27.1-foss-2018b

for i in 1 2 3 4 5
do
NAME1=`sed -n "$i p" $sample_list | awk '{print $1}'`
NAME2=`sed -n "$i p" $sample_list | awk '{print $2}'`
echo $NAME1 >$i.txt
less $NAME1 |wc -l >>$i.txt
bedtools intersect -a $NAME1 -b $NAME2 -f 0.25 -r -v |wc -l >> $i.txt
less $NAME2 |wc -l >>$i.txt
bedtools intersect -a $NAME2 -b $NAME1 -f 0.25 -r -v |wc -l >> $i.txt
done
paste *.txt>overlaps.tab
rm *.txt