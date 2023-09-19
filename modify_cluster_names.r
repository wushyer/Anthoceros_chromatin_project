setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/scale_regions")
TEcluster<- read.table("scale_regions1000.TEs.heatmap.k8.bed", header = F, sep = "\t")
#separate by cluster names
TC1<-TEcluster[grep("cluster_1", TEcluster$V13),]
TC2<-TEcluster[grep("cluster_2", TEcluster$V13),]
TC3<-TEcluster[grep("cluster_3", TEcluster$V13),]
TC4<-TEcluster[grep("cluster_4", TEcluster$V13),]
TC5<-TEcluster[grep("cluster_5", TEcluster$V13),]
TC6<-TEcluster[grep("cluster_6", TEcluster$V13),]
TC7<-TEcluster[grep("cluster_7", TEcluster$V13),]
TC8<-TEcluster[grep("cluster_8", TEcluster$V13),]
#assign new cluster names
TC1$cluster<-rep("cluster_2",5582)
TC2$cluster<-rep("cluster_8",17580)
TC3$cluster<-rep("cluster_3",20380)
TC4$cluster<-rep("cluster_1",5235)
TC5$cluster<-rep("cluster_4",5133)
TC6$cluster<-rep("cluster_5",4846)
TC7$cluster<-rep("cluster_6",10122)
TC8$cluster<-rep("cluster_7",17177)
#bind each cluster and remove old cluster names
TEcluster2<-rbind(TC4,TC1,TC3,TC5,TC6,TC7,TC8,TC2)
TEcluster3<-TEcluster2
TEcluster3$length<-TEcluster3$V3-TEcluster3$V2
TEcluster3<-TEcluster3[,c(1,2,3,4,5,6,13,14)]
TEcluster3$ID<-lapply(strsplit(TEcluster3$V4, "-"),"[[",1) #separate TE IDs
TEcluster3$ID<-unlist(TEcluster3$ID)#unlist
TEcluster3$name<-lapply(strsplit(TEcluster3$V4, "-"),"[[",2) #separate TE names
TEcluster3$name<-unlist(TEcluster3$name)#unlist
TEcluster3$family<-lapply(strsplit(TEcluster3$V4, "-"),"[[",3) #separate TE family
TEcluster3$family<-unlist(TEcluster3$family)#unlist
TEcluster3$type<-lapply(strsplit(TEcluster3$V4, "-"),"[[",4) #separate TE IDs
TEcluster3$type<-unlist(TEcluster3$type)#unlist
TEcluster3$class<-lapply(paste(TEcluster3$type,TEcluster3$family,sep = "/"),"[[",1)
TEcluster3$class<-unlist(TEcluster3$class)#unlist
TEcluster3$V1<-as.factor(TEcluster3$V1)
TEcluster3$class<-as.factor(TEcluster3$class)
TEcluster3$cluster<-as.factor(TEcluster3$cluster)
S<-TEcluster3[c(1,2,3,9,5,6,7)]
write.table(S, file = "AagrOXFv2.TE.cluster.k8.mod.all.bed",sep = "\t",row.names = F, col.names = F,quote = F)

#remove TEs in small contigs
TEcluster4<-TEcluster3[grep("AnagrOXF.S", TEcluster3$V1),]
write.table(TEcluster4, file = "AagrOXFv2.TE.cluster.k8.mod.bed",sep = "\t",row.names = F, col.names = F,quote = F)


setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS")
Genecluster<- read.table("reference_point_TSS_2000.genes.k5.bed", header = F, sep = "\t")
Genecluster$V13<-as.factor(Genecluster$V13)
summary(Genecluster$V13)
#separate by cluster names
GC1<-Genecluster[grep("cluster_1", Genecluster$V13),]
GC2<-Genecluster[grep("cluster_2", Genecluster$V13),]
GC3<-Genecluster[grep("cluster_3", Genecluster$V13),]
GC4<-Genecluster[grep("cluster_4", Genecluster$V13),]
GC5<-Genecluster[grep("cluster_5", Genecluster$V13),]
#assign new cluster names
GC1$cluster<-rep("cluster_5",length(GC1$V1))
GC2$cluster<-rep("cluster_2",length(GC2$V1))
GC3$cluster<-rep("cluster_1",length(GC3$V1))
GC4$cluster<-rep("cluster_3",length(GC4$V1))
GC5$cluster<-rep("cluster_4",length(GC5$V1))

#bind each cluster and remove old cluster names
Genecluster2<-rbind(GC3,GC2,GC4,GC5,GC1)
Genecluster2<-Genecluster2[,c(1:6,11,14)]
write.table(Genecluster2, file = "AagrOXFv2.Gene.cluster.k5.mod.all.bed",sep = "\t",row.names = F, col.names = F,quote = F)
Genecluster2$cluster<-as.factor(Genecluster2$cluster)
#remove Genes in small contigs
Genecluster3<-Genecluster2[grep("AnagrOXF.S", Genecluster2$V1),]
write.table(Genecluster3, file = "AagrOXFv2.Gene.cluster.k5.mod.bed",sep = "\t",row.names = F, col.names = F,quote = F)
