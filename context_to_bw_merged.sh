#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1G
#SBATCH --time=0-02:00:00
#SBATCH --output=logs/%x_%A_%a_%x.out
#SBATCH --error=logs/%x_%A_%a_%x.out

# post_processing
cd /groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/EMseq/merged
module load build-env/.f2021
module load  bedgraphtobigwig/385-linux-x86_64
file=(`ls *.end1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz | sed -n ${SLURM_ARRAY_TASK_ID}p`)
sample=${file%.end1_val_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz}
echo -e ". . . Sample: $sample . . .\n"

echo -e "$(date). . . creating bedgraph and bigwig files . . .\n"
for mincov in 10 
	do for j in CG CHG CHH
		do zcat $file | grep -v "mitochondrion" | grep -v "chloroplast" | awk -v threshold="$mincov" -v context="$j" 'BEGIN{OFS="\t"} $4+$5>threshold && $6==context {print $1, $2, $2+1, $4/($4+$5)}' | sort -k1,1 -k2,2n  > ${sample}_mincov${mincov}_${j}.bg
		bedGraphToBigWig ${sample}_mincov${mincov}_${j}.bg /groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2/Anthoceros_agrestis_Oxford_chrom_sizes.txt ${sample}_mincov${mincov}_${j}.bw
		done
done


