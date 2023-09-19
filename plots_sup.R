#FigureS1
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/annotations")
TE<- read.table("Aa.us.TE.bed", header = F, sep = "\t", col.names = c("chr","start","end","gene_id","score","strand"))
TE$length<-TE$end-TE$start
PCG<-read.table("AagrOXF.gene.v2.bed", header = F, sep = "\t")[,c(1,2,3,4,5,6,8)]
colnames(PCG)<-c("chr","start","end","gene_id","score","strand","length")
chr<-read.table("CHR.txt", header = F, sep = "\t", col.names = c("chr","length"))
chr$length<-chr$length/1000000
chr$TE<-c(sum(TE[grep('AnagrOXF.S1', TE$chr),]$length)/1000000,
          sum(TE[grep('AnagrOXF.S2', TE$chr),]$length)/1000000,
          sum(TE[grep('AnagrOXF.S3', TE$chr),]$length)/1000000,
          sum(TE[grep('AnagrOXF.S4', TE$chr),]$length)/1000000,
          sum(TE[grep('AnagrOXF.S5', TE$chr),]$length)/1000000,
          sum(TE[grep('AnagrOXF.S6', TE$chr),]$length)/1000000)

chr$PCG<-c(sum(PCG[grep('AnagrOXF.S1', PCG$chr),]$length)/1000000,
          sum(PCG[grep('AnagrOXF.S2', PCG$chr),]$length)/1000000,
          sum(PCG[grep('AnagrOXF.S3', PCG$chr),]$length)/1000000,
          sum(PCG[grep('AnagrOXF.S4', PCG$chr),]$length)/1000000,
          sum(PCG[grep('AnagrOXF.S5', PCG$chr),]$length)/1000000,
          sum(PCG[grep('AnagrOXF.S6', PCG$chr),]$length)/1000000)
chr$TEprop<-(chr$TE/chr$length)*100
chr$PCGprop<-(chr$PCG/chr$length)*100

chr$TEcounts<-c(length(TE[grep('AnagrOXF.S1', TE$chr),]$chr),
                length(TE[grep('AnagrOXF.S2', TE$chr),]$chr),
                length(TE[grep('AnagrOXF.S3', TE$chr),]$chr),
                length(TE[grep('AnagrOXF.S4', TE$chr),]$chr),
                length(TE[grep('AnagrOXF.S5', TE$chr),]$chr),
                length(TE[grep('AnagrOXF.S6', TE$chr),]$chr))

chr$PCGcounts<-c(length(PCG[grep('AnagrOXF.S1', PCG$chr),]$chr),
           length(PCG[grep('AnagrOXF.S2', PCG$chr),]$chr),
           length(PCG[grep('AnagrOXF.S3', PCG$chr),]$chr),
           length(PCG[grep('AnagrOXF.S4', PCG$chr),]$chr),
           length(PCG[grep('AnagrOXF.S5', PCG$chr),]$chr),
           length(PCG[grep('AnagrOXF.S6', PCG$chr),]$chr))
write.table(chr, file = "AnagrOXF_chr_summary.txt", sep = "\t",quote = F, row.names = F,col.names = T)
library(ggplot2)
#FigureS1b
ggplot(chr,aes(x=chr,y=PCGprop))+geom_bar(stat = "identity",color = "black", fill ="darkgreen")+theme_classic()+
  ylab("Proportion of PCGs")+
  xlab("Chromosome")+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14, angle = 45, hjust = 1, vjust = 1), 
        axis.title.x=element_text(size=14,face="bold",vjust = -1),
        axis.title.y=element_text(size=14,face="bold",vjust = 2))
#FigureS1c
ggplot(chr,aes(x=chr,y=TEprop))+geom_bar(stat = "identity",color = "black", fill ="tomato2")+theme_classic()+
  ylab("Proportion of PCGs")+
  xlab("Chromosome")+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14, angle = 45, hjust = 1, vjust = 1), 
        axis.title.x=element_text(size=14,face="bold",vjust = -1),
        axis.title.y=element_text(size=14,face="bold",vjust = 2))
#FigureS2
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/genmap")
data_files <- list.files("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/genmap/genmap_out")
for(i in 1:5) {                              # Head of for-loop
  assign(paste0("gcl", i),                                   # Read and store data frames
         read.table(paste0("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/genmap/genmap_out/",
                           "gcl",i,".txt")))}

for(i in 1:8) {                              # Head of for-loop
  assign(paste0("tcl", i),                                   # Read and store data frames
         read.table(paste0("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/genmap/genmap_out/",
                           "tcl",i,".txt")))
}
gcl1$cluster<-rep("cluster 1",length(gcl1$V1))
gcl2$cluster<-rep("cluster 2",length(gcl2$V1))
gcl3$cluster<-rep("cluster 3",length(gcl3$V1))
gcl4$cluster<-rep("cluster 4",length(gcl4$V1))
gcl5$cluster<-rep("cluster 5",length(gcl5$V1))
Genes<-rbind(gcl1,gcl2,gcl3,gcl4,gcl5)

tcl1$cluster<-rep("cluster 1",length(tcl1$V1))
tcl2$cluster<-rep("cluster 2",length(tcl2$V1))
tcl3$cluster<-rep("cluster 3",length(tcl3$V1))
tcl4$cluster<-rep("cluster 4",length(tcl4$V1))
tcl5$cluster<-rep("cluster 5",length(tcl5$V1))
tcl6$cluster<-rep("cluster 6",length(tcl6$V1))
tcl7$cluster<-rep("cluster 7",length(tcl7$V1))
tcl8$cluster<-rep("cluster 8",length(tcl8$V1))
TEs<-rbind(tcl1,tcl2,tcl3,tcl4,tcl5,tcl6,tcl7,tcl8)
library(ggplot2)
#FigureS2c
ggplot(Genes, aes(x=cluster,color=cluster,y=V6))+geom_boxplot()+theme_classic()+
  ylab("Mappability")+
  xlab("PCG cluster")+
  scale_color_manual(values = c("#003399","#CC3333","#003300","#009966","#666666"), labels=c("P1","P2","P3","P4","P5","P6","P7","P8"))+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.position = "none")+
  scale_x_discrete(labels=c("P1","P2","P3","P4","P5","P6","P7","P8"))

#FigureS2d
ggplot(TEs, aes(x=cluster,color=cluster,y=V6))+geom_boxplot()+theme_classic()+
  ylab("Mappability")+
  xlab("TE cluster")+
  scale_color_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33","#666666"), labels=c("T1","T2","T3","T4","T5","T6","T7","T8"))+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.position = "none")+
  scale_x_discrete(labels=c("T1","T2","T3","T4","T5","T6","T7","T8"))

#FigureS3
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/RNAseq")
exp<-read.table("exp.txt", header = T, sep = "\t")
cluster<-read.table("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS/AagrOXFv2.Gene.cluster.k5.mod.bed", header = F, sep = "\t")
#in cluster: cut -f4 GenesOverlappedbyTE25.txt|uniq>ListofGenesOverlappedByTE25.txt
OverlapTE<-read.table("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/scale_regions/ListofGenesOverlappedByTE25.txt",header = F, sep = "\t" )
colnames(cluster)<-c("chr","start","end","gene_id","score","strand","length","cluster")

exp$Ga2w_ave<-(exp$Bonn_Ga2w_rep1+exp$Bonn_Ga2w_rep2)/2
exp$Ga4w_ave<-(exp$Bonn_Ga4w_rep1+exp$Bonn_Ga4w_rep2)/2
exp$Ga8w_ave<-(exp$Bonn_Ga8w_rep1+exp$Bonn_Ga8w_rep2)/2
exp$Spsmave<-(exp$Bonn_Spsm_rep1+exp$Bonn_Spsm_rep2)/2
exp$Spmiave<-(exp$Bonn_Spmi_rep1+exp$Bonn_Spmi_rep2)/2
exp$Splaave<-(exp$Bonn_Spla_rep1+exp$Bonn_Spla_rep2)/2

df<-exp[,c(1,14:19)]
merged<-merge(cluster, df, by = "gene_id")
merged$chr<-as.factor(merged$chr)
merged$cluster<-as.factor(merged$cluster)
gene_id_overlap <- merged$gene_id %in% OverlapTE$V1
merged$overlap<-gene_id_overlap
TE<-subset(merged, merged$overlap==TRUE)
Gene<-subset(merged, merged$overlap==FALSE)

library(pheatmap)
cl1<-Gene[grep("cluster_1", Gene$cluster),]
cl1exp<-as.matrix(cl1[,9:14])
row.names(cl1exp)<-cl1$gene_id
colnames(cl1exp)<-c( "Gametophyte 2wk","Gametophyte 1mo", "Gametophyte 2mo", "Sporophyte <5mm","Sporophyte 5-10mm", "Sporophyte >10mm")
cl1exp.cluster<-kmeans(asinh(cl1exp),centers = 3)
write.table(data.frame(cl1exp.cluster$cluster), file = "CL1subcluster.txt", sep="\t", quote = F, row.names = T,col.names = T)

pheatmap(asinh(cl1exp[order(data.frame(cl1exp.cluster$cluster)),]),
         breaks = c(seq(0,10,0.1)), color = colorRampPalette(c("white","#003399"))(length(c(seq(0,10,0.1)))), 
         annotation_row = data.frame(cl1exp.cluster$cluster), cluster_rows = F, cluster_cols = F, labels_row = "")

summary(as.factor(cl1exp.cluster$cluster))

cl2<-Gene[grep("cluster_2", Gene$cluster),]
cl2exp<-as.matrix(cl2[,9:14])
row.names(cl2exp)<-cl2$gene_id
colnames(cl2exp)<-c( "Gametophyte 2wk","Gametophyte 1mo", "Gametophyte 2mo", "Sporophyte <5mm","Sporophyte 5-10mm", "Sporophyte >10mm")
cl2exp.cluster<-kmeans(asinh(cl2exp),centers = 3)
write.table(data.frame(cl2exp.cluster$cluster), file = "CL2subcluster.txt", sep="\t", quote = F, row.names = T,col.names = T)
cl2s<-read.table("CL2subcluster.txt", sep="\t",header = T)
cl2s$geneID<-row.names(cl2s)
cl2s$chr<-lapply(substr(cl2s$geneID,10,11),"[[",1)
cl2s$chr<-unlist(cl2s$chr)
cl2s$chr<-as.factor(cl2s$chr)
summary(cl2s$chr)
cl2s2<-cl2s[grep('2',cl2s$cl2exp.cluster.cluster),]
summary(cl2s2$chr)
cl2s2[grep('S6',cl2s2$chr),]
pheatmap(asinh(cl2exp[order(data.frame(cl2exp.cluster$cluster)),]),
         breaks = c(seq(0,7,0.1)), color = colorRampPalette(c("white","#CC3333"))(length(c(seq(0,7,0.1)))), 
         annotation_row = data.frame(cl2exp.cluster$cluster), cluster_rows = F, cluster_cols = F, labels_row = "")
summary(as.factor(cl2exp.cluster$cluster))

cl3<-Gene[grep("cluster_3", Gene$cluster),]
cl3exp<-as.matrix(cl3[,9:14])
row.names(cl3exp)<-cl3$gene_id
colnames(cl3exp)<-c( "Gametophyte 2wk","Gametophyte 1mo", "Gametophyte 2mo", "Sporophyte <5mm","Sporophyte 5-10mm", "Sporophyte >10mm")
cl3exp.cluster<-kmeans(asinh(cl3exp),centers = 3)
pheatmap(asinh(cl3exp[order(data.frame(cl3exp.cluster$cluster)),]),
         breaks = c(seq(0,6,0.1)), color = colorRampPalette(c("white","#003300"))(length(c(seq(0,6,0.1)))), 
         annotation_row = data.frame(cl3exp.cluster$cluster), cluster_rows = F, cluster_cols = F, labels_row = "")
summary(as.factor(cl3exp.cluster$cluster))

cl4<-Gene[grep("cluster_4", Gene$cluster),]
cl4exp<-as.matrix(cl4[,9:14])
row.names(cl4exp)<-cl4$gene_id
colnames(cl4exp)<-c( "Gametophyte 2wk","Gametophyte 1mo", "Gametophyte 2mo", "Sporophyte <5mm","Sporophyte 5-10mm", "Sporophyte >10mm")
cl4exp.cluster<-kmeans(asinh(cl4exp),centers = 3)
pheatmap(asinh(cl4exp[order(data.frame(cl4exp.cluster$cluster)),]),
         breaks = c(seq(0,10,0.1)), color = colorRampPalette(c("white","#009966"))(length(c(seq(0,10,0.1)))), 
         annotation_row = data.frame(cl4exp.cluster$cluster), cluster_rows = F, cluster_cols = F, labels_row = "")
summary(as.factor(cl4exp.cluster$cluster))

#FigureS5
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS")

closeU<-read.table("ClosestUpGenes.txt",header = F,sep="\t")
closeU$V8<-as.factor(closeU$V8)
closeU$V17<-as.factor(closeU$V17)
closeU<-closeU[grep('cluster',closeU$V17),]
closeU<-closeU[grep('cluster',closeU$V8),]
closeUGenes<-closeU[grep('1',closeU$V9),]
closeUGenes1<-closeUGenes[grep('cluster_5',closeUGenes$V8,invert = T),]
closeUGenes1<-closeUGenes1[grep('cluster_5',closeUGenes1$V17,invert = T),]
closeUTEs<-closeU[grep('2',closeU$V9),]
closeUTEs1<-closeUTEs[grep('cluster_5',closeUTEs$V8,invert = T),]
closeUTEs1<-closeUTEs1[grep('cluster_8',closeUTEs1$V17,invert = T),]

closeD<-read.table("ClosestDownGenes.txt",header = F,sep="\t")
closeD$V8<-as.factor(closeD$V8)
closeD$V17<-as.factor(closeD$V17)
closeD<-closeD[grep('cluster',closeD$V17),]
closeD<-closeD[grep('cluster',closeD$V8),]
closeDGenes<-closeD[grep('1',closeD$V9),]
closeDGenes1<-closeDGenes[grep('cluster_5',closeDGenes$V8,invert = T),]
closeDGenes1<-closeDGenes1[grep('cluster_5',closeDGenes1$V17,invert = T),]
closeDTEs<-closeD[grep('2',closeD$V9),]
closeDTEs1<-closeDTEs[grep('cluster_5',closeDTEs$V8,invert = T),]
closeDTEs1<-closeDTEs1[grep('cluster_8',closeDTEs1$V17,invert = T),]
#free_y histogram
ggplot(closeUGenes1, aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0))+
  facet_grid(rows='V8',scales = "free_y")+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeUTEs1, aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0))+
  facet_grid(rows='V8',scales = "free_y")+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))

ggplot(closeDGenes1, aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0,5000))+
  facet_grid(rows='V8',scales = "free_y")+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeDTEs1, aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0,5000))+
  facet_grid(rows='V8',scales = "free_y")+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))

closeTU<-read.table("ClosestUpTEs.txt", sep = "\t", header = F)
closeTU$V8<-as.factor(closeTU$V8)
closeTU$V17<-as.factor(closeTU$V17)
closeTU<-closeTU[grep('cluster',closeTU$V17),]
closeTU<-closeTU[grep('cluster',closeTU$V8),]
closeTUGenes<-closeTU[grep('1',closeTU$V9),]
closeTUGenes1<-closeTUGenes[grep('cluster_8',closeTUGenes$V8,invert = T),]
closeTUGenes1<-closeTUGenes1[grep('cluster_5',closeTUGenes1$V17,invert = T),]
closeTUTEs<-closeTU[grep('2',closeTU$V9),]
closeTUTEs1<-closeTUTEs[grep('cluster_8',closeTUTEs$V8,invert = T),]
closeTUTEs1<-closeTUTEs1[grep('cluster_8',closeTUTEs1$V17,invert = T),]

closeTD<-read.table("ClosestDownTEs.txt", sep = "\t", header = F)
closeTD$V8<-as.factor(closeTD$V8)
closeTD$V17<-as.factor(closeTD$V17)
closeTD<-closeTD[grep('cluster',closeTD$V17),]
closeTD<-closeTD[grep('cluster',closeTD$V8),]
closeTDGenes<-closeTD[grep('1',closeTD$V9),]
closeTDGenes1<-closeTDGenes[grep('cluster_8',closeTDGenes$V8,invert = T),]
closeTDGenes1<-closeTDGenes1[grep('cluster_5',closeTDGenes1$V17,invert = T),]
closeTDTEs<-closeTD[grep('2',closeTD$V9),]
closeTDTEs1<-closeTDTEs[grep('cluster_8',closeTDTEs$V8,invert = T),]
closeTDTEs1<-closeTDTEs1[grep('cluster_8',closeTDTEs1$V17,invert = T),]
ggplot(closeTUGenes1, aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0))+
  facet_grid(rows='V8',scales = "free_y")+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeTUTEs1, aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0))+
  facet_grid(rows='V8',scales = "free_y")+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))

ggplot(closeTDGenes1, aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0, 5000))+
  facet_grid(rows='V8',scales = "free_y")+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeTDTEs1, aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0, 5000))+
  facet_grid(rows='V8',scales = "free_y")+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))

