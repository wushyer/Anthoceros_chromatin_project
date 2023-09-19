bedtools intersect -b AagrOXFv2.Gene.cluster.k5.mod.sorted.bed  -a AagrOXFv2.TE.cluster.k8.mod.soreted.bed -f 0.25 -wo >TEsOverlappedbyGenes25.txt
bedtools intersect -a AagrOXFv2.Gene.cluster.k5.mod.sorted.bed  -b AagrOXFv2.TE.cluster.k8.mod.soreted.bed -f 0.25 -wo >GenesOverlappedbyTE25.txt
