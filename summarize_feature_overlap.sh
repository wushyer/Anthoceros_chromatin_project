#!/bin/bash

# === begin SBATCH directives ===
#SBATCH --partition=c
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=1-5
#SBATCH --output=logs/%x_%A_%a.txt
# === end SBATCH directives ===
#Summarize ratio of overlapping bp of peak and features
#Created 8.1.2018 by Sean Montgomery
#modified 07.09.2022 by Tetsuya Hisanaga

#sort gff using grep'chromosome coordinate' input.gff | sort -k1,1V -k4,4n

sample_list=/groups/berger/lab/tetsuya.hisanaga/input_files/20220907.txt 
work_folder=/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam
#sample name
F=`sed -n "${SLURM_ARRAY_TASK_ID} p" $sample_list | awk '{print $1}'`

cd $work_folder
module load bedtools/2.27.1-foss-2018b

#Sample usage:
# for f in *bedgraph; do /volumes/berger/user/sean.montgomery/Documents/Scripts/chipseq/summarize_feature_overlap.sh $f; done
 #1. Sample bedgraph file from macs2 broadpeak
 file=macs2/${F}_peaks.broadPeak
 cut -f 1-4 $file > ${F}.bed
 file1=${F}.bed
 #2 gff3 file of genomic features
 file2=/groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2/Anthoceros_agrestis_Oxford_gene_annotations.sorted.gff
 file3=/groups/berger/lab/tetsuya.hisanaga/ref_files/Anthoceros/AagrOXFv2/Aa.us.fa.mod.EDTA.TEanno.gff3
 # Summarize feature overlap
 # feature=$2 #CDS, five_prime_UTR, three_prime_UTR, exon, mRNA, intron
 # coverage_length=`bedtools intersect -a $file1 -b $file2 -wo | grep $feature | awk '{sum+=$14} END {print sum}'`
 # feature_length=`grep $feature $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 ip_length=`awk '{a=$3-$2;print $0,a;}' $file1 | awk '{sum+=$5} END {print sum}'`

 cov_exon=`bedtools intersect -a $file1 -b $file2 -wo | grep -E '\bexon\b' | grep '\.t1;' | awk '{sum+=$14} END {print sum}'`
 cov_intron=`bedtools intersect -a $file1 -b $file2 -wo | grep -E '\bintron\b' | grep '\.t1;' | awk '{sum+=$14} END {print sum}'`
 cov_gene=`bedtools intersect -a $file1 -b $file2 -wo | grep -E '\bgene\b' | awk '{sum+=$14} END {print sum}'`
 cov_CACTA=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bCACTA_TIR_transposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_hAT=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bhAT_TIR_transposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_PIF=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bPIF_Harbinger_TIR_transposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_Tc1=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bTc1_Mariner_TIR_transposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_Mutator=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bMutator_TIR_transposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_helitron=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bhelitron\b' | awk '{sum+=$14} END {print sum}'`
 cov_Copia=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bCopia_LTR_retrotransposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_Gypsy=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bGypsy_LTR_retrotransposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_LTR=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bLTR_retrotransposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_LINE=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bLINE_element\b' | awk '{sum+=$14} END {print sum}'`
 cov_penelope=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\bPenelope_retrotransposon\b' | awk '{sum+=$14} END {print sum}'`
 cov_repeat=`bedtools intersect -a $file1 -b $file3 -wo | grep -E '\brepeat_region\b' | awk '{sum+=$14} END {print sum}'`
 
 exon=`grep -E '\bexon\b' $file2 |grep '\.t1;' | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 intron=`grep -E '\bintron\b' $file2 | grep '\.t1;' | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 gene=`grep -E '\bgene\b'  $file2 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 CACTA=`grep -E '\bCACTA_TIR_transposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 hAT=`grep -E '\bhAT_TIR_transposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 PIF=`grep -E '\bPIF_Harbinger_TIR_transposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 Tc1=`grep -E '\bTc1_Mariner_TIR_transposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 Mutator=`grep -E '\bMutator_TIR_transposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 helitron=`grep -E '\bhelitron\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 Copia=`grep -E '\bCopia_LTR_retrotransposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 Gypsy=`grep -E '\bGypsy_LTR_retrotransposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 LTR=`grep -E '\bLTR_retrotransposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 LINE=`grep -E '\bLINE_element\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 penelope=`grep -E '\bPenelope_retrotransposon\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 repeat=`grep -E '\brepeat_region\b' $file3 | awk '{a=$5-$4;print $0,a;}' | awk '{sum+=$10} END {print sum}'`
 
 UTR=`echo "scale=3; $gene - $exon - $intron" | bc -l`
 cov_UTR=`echo "scale=3; $cov_gene - $cov_exon - $cov_intron" | bc -l`
 
 echo "Total length of" $file1 ":" $ip_length

 echo "Total coverage length for UTR in" $file1 ":" $cov_UTR
 UTRer=`echo "scale=3; $cov_UTR / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for UTR:" $UTR
 UTRer2=`echo "scale=3; $cov_UTR / $UTR" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for exon in" $file1 ":" $cov_exon
 exoner=`echo "scale=3; $cov_exon / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for exon:" $exon
 exoner2=`echo "scale=3; $cov_exon / $exon" | bc -l` #Portion of feature that is sample
 echo ""
  echo "Total coverage length for intron in" $file1 ":" $cov_intron
 introner=`echo "scale=3; $cov_intron / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for intron:" $intron
 introner2=`echo "scale=3; $cov_intron / $intron" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for repeats in" $file1 ":" $cov_repeat
 repeater=`echo "scale=3; $cov_repeat / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for repeats:" $repeat
 repeater2=`echo "scale=3; $cov_repeat / $repeat" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for CACTA in" $file1 ":" $cov_CACTA
 cactaer=`echo "scale=3; $cov_CACTA / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for CACTA:" $CACTA
 cactaer2=`echo "scale=3; $cov_CACTA / $CACTA" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for hAT in" $file1 ":" $cov_hAT
 hater=`echo "scale=3; $cov_hAT / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for hAT:" $hAT
 hater2=`echo "scale=3; $cov_hAT / $hAT" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for PIF in" $file1 ":" $cov_PIF
 pifer=`echo "scale=3; $cov_PIF / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for PIF:" $PIF
 pifer2=`echo "scale=3; $cov_PIF / $PIF" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for Tc1 in" $file1 ":" $cov_Tc1
 tc1er=`echo "scale=3; $cov_Tc1 / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for Tc1:" $Tc1
 tc1er2=`echo "scale=3; $cov_Tc1 / $Tc1" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for Mutator in" $file1 ":" $cov_Mutator
 mutatorer=`echo "scale=3; $cov_Mutator / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for Mutator:" $Mutator
 mutatorer2=`echo "scale=3; $cov_Mutator / $Mutator" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for helitron in" $file1 ":" $cov_helitron
 helier=`echo "scale=3; $cov_helitron / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for helitron:" $helitron
 helier2=`echo "scale=3; $cov_helitron / $helitron" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for Copia in" $file1 ":" $cov_Copia
 Copiaer=`echo "scale=3; $cov_Copia / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for Copia:" $Copia
 Copiaer2=`echo "scale=3; $cov_Copia / $Copia" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for Gypsy in" $file1 ":" $cov_Gypsy
 Gypsyer=`echo "scale=3; $cov_Gypsy / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for Gypsy:" $Gypsy
 Gypsyer2=`echo "scale=3; $cov_Gypsy / $Gypsy" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for LTR in" $file1 ":" $cov_LTR
 LTRer=`echo "scale=3; $cov_LTR / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for LTR:" $LTR
 LTRer2=`echo "scale=3; $cov_LTR / $LTR" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for LINE in" $file1 ":" $cov_LINE
 LINEer=`echo "scale=3; $cov_LINE / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for LINE:" $LINE
 LINEer2=`echo "scale=3; $cov_LINE / $LINE" | bc -l` #Portion of feature that is sample
 echo ""
 echo "Total coverage length for penelope in" $file1 ":" $cov_penelope
 penelopeer=`echo "scale=3; $cov_penelope / $ip_length" | bc -l` #Portion of sample that is feature
 echo "Total feature length for penelope:" $penelope
 penelopeer2=`echo "scale=3; $cov_penelope / $penelope" | bc -l` #Portion of feature that is sample
 echo ""
  NAME=`basename $file1 .bed`
 # echo -e "Sample \t 5primeUTR \t CDS \t 3primeUTR \t intron \t repeats \t gene \t rRNA \t gene \t simple \t satellite \t unknown \t repeat" > sample_overlap_summary.txt
 echo -e "$NAME \t $UTRer \t $exoner \t  $introner \t $repeater \t $cactaer \t $hater \t $pifer \t $tc1er \t $mutatorer \t $helier \t $Copiaer \t $Gypsyer \t $LTRer \t $LINEer \t $penelopeer"  >> sample_overlap_summary.txt
 # echo -e "Sample \t 5primeUTR \t CDS \t 3primeUTR \t intron \t repeats \t gene \t rRNA \t gene \t simple \t satellite \t unknown \t repeat" > feature_overlap_summary.txt
 echo -e "$NAME \t $UTRer2 \t $exoner2 \t $introner2 \t $repeater2 \t $cactaer2 \t $hater2 \t $pifer2 \t $tc1er2\t $mutatorer2 \t $helier2 \t $Copiaer2 \t $Gypsyer2 \t $LTRer2 \t $LINEer2 \t $penelopeer2" >> feature_overlap_summary.txt
# #By clusters
# #1. Determine genes in clusters
# # grep 'cluster 1' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster1.txt
# # grep 'cluster 2' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster2.txt
# # grep 'cluster 3' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster3.txt
# # grep 'cluster 4' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster4.txt
# # grep 'cluster 5' selected.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.TSS.profiledata.k5.heatmap.cluster5.txt
# #1. Sample bedgraph file from homer
# file1=$1
# # Summarize feature overlap
# ip_length=`awk '{a=$3-$2;print $0,a;}' $file1 | awk '{sum+=$5} END {print sum}'`
# for i in {1..5}; do
# 	file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.cluster$i.txt
# 	cov_gene=`bedtools intersect -a $file1 -b $file2 -wo | awk '{sum+=$13} END {print sum}'`
# 	gene=`awk '{a=$3-$2;print $0,a;}' $file2 | awk '{sum+=$9} END {print sum}'`
# 	echo "Total coverage length for genes in" $file1 ":" $cov_gene
# 	echo "Total feature length for genes:" $gene
# 	echo ""
# 	eval "samplecluster$i=`echo "scale=3; $cov_gene / $ip_length" | bc -l`"
# 	eval "featurecluster$i=`echo "scale=3; $cov_gene / $gene" | bc -l`"
# done
# echo "Total length of" $file1 ":" $ip_length
# NAME=`basename $file1 _Tak1v4.regions.narrow.bedgraph`
# # echo -e "Sample \t Cluster_1 \t Cluster_2 \t Cluster_3 \t Cluster_4 \t Cluster_5" > sample_overlap_summary.clusters.txt
# echo -e "$NAME \t $samplecluster1 \t $samplecluster2 \t $samplecluster3 \t $samplecluster4 \t $samplecluster5" >> sample_overlap_summary.clusters.txt
# # echo -e "Sample \t Cluster_1 \t Cluster_2 \t Cluster_3 \t Cluster_4 \t Cluster_5" > feature_overlap_summary.clusters.txt
# echo -e "$NAME \t $featurecluster1 \t $featurecluster2 \t $featurecluster3 \t $featurecluster4 \t $featurecluster5" >> feature_overlap_summary.clusters.txt
# ##By quintile
# ##1. Determine quintiles
# # quintiles <- quantile(gametop$mean, probs = seq(0,1,0.2))
# # bedlist <- read.table("/volumes/berger/lab/NGS_annotations/backup_berger_common/Tak1v4/annotations/Tak1v4.gene.bed")
# # bedlist$TPM <- gametop$mean
# # write.table(subset(bedlist,bedlist$TPM<quintiles[[2]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile1.txt",sep='\t',col.names = F,row.names = F)
# # write.table(subset(bedlist,bedlist$TPM>=quintiles[[5]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile5.txt",sep='\t',col.names = F,row.names = F)
# # write.table(subset(bedlist,bedlist$TPM>=quintiles[[4]] & bedlist$TPM<quintiles[[5]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile4.txt",sep='\t',col.names = F,row.names = F)
# # write.table(subset(bedlist,bedlist$TPM>=quintiles[[3]] & bedlist$TPM<quintiles[[4]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile3.txt",sep='\t',col.names = F,row.names = F)
# # write.table(subset(bedlist,bedlist$TPM>=quintiles[[2]] & bedlist$TPM<quintiles[[3]]), file="/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile2.txt",sep='\t',col.names = F,row.names = F)
# ##2. Sample bedgraph file from homer
# file1=$1
# # Summarize feature overlap
# ip_length=`awk '{a=$3-$2;print $0,a;}' $file1 | awk '{sum+=$5} END {print sum}'`
# for i in {1..5}; do
# 	file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/expression/homer/narrow/gametop.expr.quintile$i.txt
# 	sed -i -e 's/"//g' $file2
# 	cov_gene=`bedtools intersect -a $file1 -b $file2 -wo | awk '{sum+=$11} END {print sum}'`
# 	gene=`awk '{a=$3-$2;print $0,a;}' $file2 | awk '{sum+=$7} END {print sum}'`
# 	echo "Total coverage length for genes in" $file1 ":" $cov_gene
# 	echo "Total feature length for genes:" $gene
# 	echo ""
# 	eval "samplecluster$i=`echo "scale=3; $cov_gene / $ip_length" | bc -l`"
# 	eval "featurecluster$i=`echo "scale=3; $cov_gene / $gene" | bc -l`"
# done
# echo "Total length of" $file1 ":" $ip_length
# NAME=`basename $file1 _Tak1v4.regions.narrow.bedgraph`
# # echo -e "Sample \t Quintile_1 \t Quintile_2 \t Quintile_3 \t Quintile_4 \t Quintile_5" > sample_overlap_summary.quintiles.txt
# echo -e "$NAME \t $samplecluster1 \t $samplecluster2 \t $samplecluster3 \t $samplecluster4 \t $samplecluster5" >> sample_overlap_summary.quintiles.txt
# # echo -e "Sample \t Quintile_1 \t Quintile_2 \t Quintile_3 \t Quintile_4 \t Quintile_5" > feature_overlap_summary.quintiles.txt
# echo -e "$NAME \t $featurecluster1 \t $featurecluster2 \t $featurecluster3 \t $featurecluster4 \t $featurecluster5" >> feature_overlap_summary.quintiles.txt
# ##4. Summarize repeat cluster overlap
# echo -e "Cluster \t Quintile1 \t Quintile2 \t Quintile3 \t Quintile4 \t Quintile5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/sample_overlap_summary.atac.quintiles.cluster.repeat.txt
# echo -e "Cluster \t Quintile1 \t Quintile2 \t Quintile3 \t Quintile4 \t Quintile5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/feature_overlap_summary.atac.quintiles.cluster.repeat.txt
# for i in {1..5}; do
# 	file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.cluster$i.bed.txt
# 	ip_length=`awk '{a=$3-$2;print $0,a;}' $file1 | awk '{sum+=$9} END {print sum}'`
# 	for j in {1..5}; do
# 		file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/atac.quintile$j.txt
# 		cov_gene=`bedtools intersect -a $file1 -b $file2 -wo | awk '{sum+=$19} END {print sum}'`
# 		gene=`awk '{a=$3-$2;print $0,a;}' $file2 | awk '{sum+=$11} END {print sum}'`
# 		echo "Total coverage length for genes in" $file1 ":" $cov_gene
# 		echo "Total feature length for genes:" $gene
# 		echo ""
# 		eval "samplecluster$j=`echo "scale=3; $cov_gene / $ip_length" | bc -l`"
# 		eval "featurecluster$j=`echo "scale=3; $cov_gene / $gene" | bc -l`"
# 	done
# 	NAME=`echo "Repeat"$i`
# 	echo -e "$NAME \t $samplecluster1 \t $samplecluster2 \t $samplecluster3 \t $samplecluster4 \t $samplecluster5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/sample_overlap_summary.atac.quintiles.cluster.repeat.txt
# 	echo -e "$NAME \t $featurecluster1 \t $featurecluster2 \t $featurecluster3 \t $featurecluster4 \t $featurecluster5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/atac/feature_overlap_summary.atac.quintiles.cluster.repeat.txt
# done
# ##Overlap between peaks of different marks
# ## Ratio of mark (rowname) in (colname) (overlap divided by total length of rowname)
# # Summarize feature overlap
# ls *bedgraph | sed 's/_Tak1v4.regions.narrow.bedgraph//g' | tr '\n' '\t' > sample_overlap_summary.peaks.txt
# echo "" >> sample_overlap_summary.peaks.txt
# for f in *bedgraph; do
# 	NAME=`basename $f _Tak1v4.regions.narrow.bedgraph`
# 	eval "output1=`echo $NAME`"
# 	for file in *bedgraph; do
# 		f_ip_length=`awk '{a=$3-$2;print $0,a;}' $f | awk '{sum+=$5} END {print sum}'`
# 		cov_gene=`bedtools intersect -a $f -b $file -wo | awk '{sum+=$9} END {print sum}'`
# 		echo "Total coverage length for " $f " in " $file ":" $cov_gene
# 		echo ""
# 		output1+="boop"
# 		eval "output1+=`echo "scale=3; $cov_gene / $f_ip_length" | bc -l`"
# 	done
# 	echo $output1 | sed 's/boop/	/g' >> sample_overlap_summary.peaks.txt
# done
# sed -i 's/1.000/0.000/g' sample_overlap_summary.peaks.txt
# # ##Distance between features in certain clusters
# # ##1. Determine repeat clusters
# # # grep 'cluster 1' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster1.txt
# # # grep 'cluster 2' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster2.txt
# # # grep 'cluster 3' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster3.txt
# # # grep 'cluster 4' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster4.txt
# # # grep 'cluster 5' selected.repeat.TSS.profiledata.k5.heatmap.txt | grep -v '#' | sort -k1,1 -k2,2n > selected.repeat.TSS.profiledata.k5.heatmap.cluster5.txt
# # # Find average distance of nearest gene in specified cluster to repeat in specified cluster
# # echo -e "Cluster \t Gene_1 \t Gene_2 \t Gene_3 \t Gene_4 \t Gene_5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.txt
# # echo -e "Cluster \t Gene_1 \t Gene_2 \t Gene_3 \t Gene_4 \t Gene_5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.txt
# # for i in {1..5}; do
# # 	for j in {1..5}; do
# # 		file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.cluster$i.bed.txt
# # 		file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.cluster$j.bed.txt
# # 		eval "distance$j=`bedtools closest -a $file1 -b $file2 -d | awk '{ total += $17 } END { print total/NR }'`"
# # 		eval "distancegene$j=`bedtools closest -a $file2 -b $file1 -d | awk '{ total += $17 } END { print total/NR }'`"
# # 	done
# # 	echo -e "Repeat_$i \t $distance1 \t $distance2 \t $distance3 \t $distance4 \t $distance5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.txt
# # 	echo -e "Repeat_$i \t $distancegene1 \t $distancegene2 \t $distancegene3 \t $distancegene4 \t $distancegene5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.txt
# # done
# echo -e "Repeat1.Gene1 \t Repeat1.Gene2 \t Repeat1.Gene3 \t Repeat1.Gene4 \t Repeat1.Gene5 \t Repeat2.Gene1 \t Repeat2.Gene2 \t Repeat2.Gene3 \t Repeat2.Gene4 \t Repeat2.Gene5 \t Repeat3.Gene1 \t Repeat3.Gene2 \t Repeat3.Gene3 \t Repeat3.Gene4 \t Repeat3.Gene5 \t Repeat4.Gene1 \t Repeat4.Gene2 \t Repeat4.Gene3 \t Repeat4.Gene4 \t Repeat4.Gene5 \t Repeat5.Gene1 \t Repeat5.Gene2 \t Repeat5.Gene3 \t Repeat5.Gene4 \t Repeat5.Gene5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.all.txt
# echo -e "Gene1.Repeat1 \t Gene1.Repeat2 \t Gene1.Repeat3 \t Gene1.Repeat4 \t Gene1.Repeat5 \t Gene2.Repeat1 \t Gene2.Repeat2 \t Gene2.Repeat3 \t Gene2.Repeat4 \t Gene2.Repeat5 \t Gene3.Repeat1 \t Gene3.Repeat2 \t Gene3.Repeat3 \t Gene3.Repeat4 \t Gene3.Repeat5 \t Gene4.Repeat1 \t Gene4.Repeat2 \t Gene4.Repeat3 \t Gene4.Repeat4 \t Gene4.Repeat5 \t Gene5.Repeat1 \t Gene5.Repeat2 \t Gene5.Repeat3 \t Gene5.Repeat4 \t Gene5.Repeat5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.all.txt
# for i in {1..5}; do
# 	for j in {1..5}; do
# 		file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.cluster$i.bed.txt
# 		file2=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.TSS.profiledata.k5.heatmap.cluster$j.bed.txt
# 		bedtools closest -a $file1 -b $file2 -d | cut -f15 > repeat$i.gene$j
# 		bedtools closest -a $file2 -b $file1 -d | cut -f15 > gene$j.repeat$i
# 	done
# done 
# paste repeat1.gene1 repeat1.gene2 repeat1.gene3 repeat1.gene4 repeat1.gene5 repeat2.gene1 repeat2.gene2 repeat2.gene3 repeat2.gene4 repeat2.gene5 repeat3.gene1 repeat3.gene2 repeat3.gene3 repeat3.gene4 repeat3.gene5 repeat4.gene1 repeat4.gene2 repeat4.gene3 repeat4.gene4 repeat4.gene5 repeat5.gene1 repeat5.gene2 repeat5.gene3 repeat5.gene4 repeat5.gene5 >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_to_gene_distance.clusters.all.txt
# paste gene1.repeat1 gene1.repeat2 gene1.repeat3 gene1.repeat4 gene1.repeat5 gene2.repeat1 gene2.repeat2 gene2.repeat3 gene2.repeat4 gene2.repeat5 gene3.repeat1 gene3.repeat2 gene3.repeat3 gene3.repeat4 gene3.repeat5 gene4.repeat1 gene4.repeat2 gene4.repeat3 gene4.repeat4 gene4.repeat5 gene5.repeat1 gene5.repeat2 gene5.repeat3 gene5.repeat4 gene5.repeat5 >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/gene_to_repeat_distance.clusters.all.txt
# ##Number of repeat class in a cluster
# ##1. Determine repeat clusters
# echo -e "Cluster \t Unknown \t LTR \t DNA \t LINE \t Helitron" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.txt
# echo -e "Cluster \t Unknown \t LTR \t DNA \t LINE \t Helitron" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.txt
# file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.ids.txt
# unknown_total=`grep "Unknown" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# ltr_total=`grep "LTR" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# dna_total=`grep "DNA" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# line_total=`grep "LINE" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# helitron_total=`grep "Helitron" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# for i in {1..5}; do
# 	cluster_total=`grep "cluster $i" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	unknown_length=`grep "Unknown" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	ltr_length=`grep "LTR" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	dna_length=`grep "DNA" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	line_length=`grep "LINE" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	helitron_length=`grep "Helitron" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	eval "unknown_ratio=`echo "scale=3; $unknown_length / $unknown_total" | bc -l`"
# 	eval "ltr_ratio=`echo "scale=3; $ltr_length / $ltr_total" | bc -l`"
# 	eval "dna_ratio=`echo "scale=3; $dna_length / $dna_total" | bc -l`"
# 	eval "line_ratio=`echo "scale=3; $line_length / $line_total" | bc -l`"
# 	eval "helitron_ratio=`echo "scale=3; $helitron_length / $helitron_total" | bc -l`"
# 	echo -e "Repeat_$i \t $unknown_ratio \t $ltr_ratio \t $dna_ratio \t $line_ratio \t $helitron_ratio" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.txt
# 	eval "unknown_ratio=`echo "scale=3; $unknown_length / $cluster_total" | bc -l`"
# 	eval "ltr_ratio=`echo "scale=3; $ltr_length / $cluster_total" | bc -l`"
# 	eval "dna_ratio=`echo "scale=3; $dna_length / $cluster_total" | bc -l`"
# 	eval "line_ratio=`echo "scale=3; $line_length / $cluster_total" | bc -l`"
# 	eval "helitron_ratio=`echo "scale=3; $helitron_length / $cluster_total" | bc -l`"
# 	echo -e "Repeat_$i \t $unknown_ratio \t $ltr_ratio \t $dna_ratio \t $line_ratio \t $helitron_ratio" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.txt
# done
# ##TE superfamilies by repeat cluster overlap
# echo -e "Repeat \t Cluster1 \t Cluster2 \t Cluster3 \t Cluster4 \t Cluster5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.all.txt
# echo -e "Repeat \t Cluster1 \t Cluster2 \t Cluster3 \t Cluster4 \t Cluster5" > /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.all.txt
# file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.ids.txt
# while read p; do
# 	file1=/volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/selected.repeat.TSS.profiledata.k5.heatmap.ids.txt
# 	repeat_total=`grep "$p" $file1 | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 	for i in {1..5}; do
# 		cluster_total=`grep "cluster $i" $file1 | grep -v "Unknown" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 		repeat_length=`grep "$p" $file1 | grep "cluster $i" | awk '{a=$3-$2;print $0,a;}' | wc -l`
# 		eval "ratioperclass$i=`echo "scale=3; $repeat_length / $repeat_total" | bc -l`"
# 		eval "ratiopercluster$i=`echo "scale=3; $repeat_length / $cluster_total" | bc -l`"
# 	done
# 	echo -e "$p \t $ratioperclass1 \t $ratioperclass2 \t $ratioperclass3 \t $ratioperclass4 \t $ratioperclass5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratioperclass.all.txt
# 	echo -e "$p \t $ratiopercluster1 \t $ratiopercluster2 \t $ratiopercluster3 \t $ratiopercluster4 \t $ratiopercluster5" >> /volumes/berger/user/sean.montgomery/Documents/cutrun/Tak1v4/merged/profiles/repeat_classes.clusters.ratiopercluster.all.txt	
# done </volumes/berger/lab/Marchantia_v4.4/repeat.superfamilies.txt 
