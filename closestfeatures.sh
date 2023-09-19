bedtools intersect -a AagrOXFv2.Gene.cluster.k5.mod.sorted.bed -b AagrOXFv2.TE.cluster.k8.mod.soreted.bed -v >Genes.k5.wo.TEs.bed
bedtools intersect -b AagrOXFv2.Gene.cluster.k5.mod.sorted.bed -a AagrOXFv2.TE.cluster.k8.mod.soreted.bed -v >TEs.k8.wo.Genes.bed

cut -f1,2,3,4,5,6 TEs.k8.wo.Genes.bed >tmp
cut -f8 TEs.k8.wo.Genes.bed >tmp2
cut -f7 TEs.k8.wo.Genes.bed >tmp3
paste tmp tmp2 tmp3>TEs.k8.wo.Genes.bed
bedtools closest -a Genes.k5.wo.TEs.bed -b Genes.k5.wo.TEs.bed TEs.k8.wo.Genes.bed -D a -io -id -mdb all -t first >ClosestUpGenes.txt
bedtools closest -a TEs.k8.wo.Genes.bed -b Genes.k5.wo.TEs.bed TEs.k8.wo.Genes.bed -D a -io -id -mdb all -t first >ClosestUpTEs.txt

bedtools closest -a Genes.k5.wo.TEs.bed -b Genes.k5.wo.TEs.bed TEs.k8.wo.Genes.bed -D a -io -iu -mdb all -t first >ClosestDownGenes.txt
bedtools closest -a TEs.k8.wo.Genes.bed -b Genes.k5.wo.TEs.bed TEs.k8.wo.Genes.bed -D a -io -iu -mdb all -t first >ClosestDownTEs.txt