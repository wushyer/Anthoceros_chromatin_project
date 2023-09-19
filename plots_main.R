#Figure1a overlap between chromatin marks and genomic features
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam")
OL1<-read.table("sample_overlap_summary.txt", header = F, sep = "\t")
colnames(OL1)<-c("Mark","UTR","exon","intron","other repeat","CACTA_TIR_transposon","hAT_TIR_transposon","PIF_Harbinger_TIR_transposon","Tc1_Mariner_TIR_transposon","Mutator_TIR_transposon","helitron","Copia_LTR_retrotransposon","Gypsy_LTR_retrotransposon","other_LTR_retrotransposon","LINE_element","Penelope_retrotransposon")
library(reshape2)
OL2<-melt(OL1[,c(1:6,10:14)],id.vars = "Mark", variable.name = "feature")
library(tidyverse)
levels(OL2$Mark)
g<-OL2 %>%mutate(Mark = fct_relevel(Mark, "H3K9me1 ", "H3K27me1 ","H3K4me3 ","H3K36me3 ","H3K27me3 "))
library(ggplot2)
ggplot(g[grep('UTR',g$feature,invert = T),],aes(x=Mark,y=value,fill=feature))+
  geom_bar(stat = "identity",color = "white")+
  theme_classic()+
  scale_fill_manual(values = c("palegreen3","chartreuse2","#330000","sienna4","sienna1","salmon","orange4","orange3","orange2","red4"))+
  ylab("Propotion of mark overlapping genomic feature")+
  xlab("Mark")+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14, vjust = 1, hjust = 1,angle = 45), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

#Figure2c PCG expression per chromosome
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/RNAseq")
exp<-read.table("exp.txt", header = T, sep = "\t")
cluster<-read.table("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS/AagrOXFv2.Gene.cluster.k5.mod.bed", header = F, sep = "\t")
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

library(ggplot2)
#plot geneexpression per chromosome
ggplot(na.omit(merged), aes(x=chr,y=asinh(Ga4w_ave)))+geom_violin(fill="white")+theme_classic()+
  ylab("asinh(TPM)")+
  xlab("")+
  theme(axis.text.y=element_text(angle=90,size=14),
        axis.text.x=element_text(angle=90,size=14), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2))+
  stat_summary(fun=median, geom="point", size=2, color="red")

#Figure3a
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/macs2")
peaks<-as.matrix(read.table("overlaps.tab", header = T, fill = T, sep = "\t", check.names = F))

exp1<-read.table("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/RNAseq/Bonn_Ga4w_rep1/Bonn_Ga4w_rep1.genes.results", header = T, sep = "\t")
exp2<-read.table("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/RNAseq/Bonn_Ga4w_rep2/Bonn_Ga4w_rep2.genes.results", header = T, sep = "\t")
exp<-exp1[,c(1,5,6,7)]
exp$expected_count<-(exp1$expected_count+exp2$expected_count)/2
exp$TPM<-(exp1$TPM+exp2$TPM)/2
exp$FPKM<-(exp1$FPKM+exp2$FPKM)/2

for (i in 1:length(colnames(peaks))) {
  exp[,paste0(colnames(peaks)[i],"_peak")] <- "No"
  exp[,paste0(colnames(peaks)[i],"_peak")][exp$gene_id %in% peaks[,i]] <- "Yes"
}
K9me1<-subset(exp,exp$H3K9me1_peak=="Yes")
K9me1$peak<-rep("H3K9me1",length(rownames(K9me1)))
K27me1<-subset(exp,exp$H3K27me1_peak=="Yes")
K27me1$peak<-rep("H3K27me1",length(rownames(K27me1)))
K4me3<-subset(exp,exp$H3K4me3_peak=="Yes")
K4me3$peak<-rep("H3K4me3",length(rownames(K4me3)))
K36me3<-subset(exp,exp$H3K36me3_peak=="Yes")
K36me3$peak<-rep("H3K36me3",length(rownames(K36me3)))
K27me3<-subset(exp,exp$H3K27me3_peak=="Yes")
K27me3$peak<-rep("H3K27me3",length(rownames(K27me3)))

df<-rbind(K9me1,K27me1,K4me3,K36me3,K27me3)
df<-df[,c(1,3,10)]
library(ggplot2)
library(tidyverse)
g<-df %>%mutate(peak = fct_relevel(peak, "H3K9me1", "H3K27me1", "H3K4me3", "H3K36me3", "H3K27me3"))%>%
  ggplot(aes(x=peak,y=asinh(TPM),color=peak))+geom_violin()
g+xlab("Mark")+ylab("asinh(TPM)")+
  theme_classic()+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(angle = 45,size=14,hjust = 1,vjust = 1), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.position = "none")+
  scale_color_manual(values = c("#CC6600","#CC9900","#339900","#009933","#0066CC"))+
  stat_summary(fun.y=median, geom="point", size=2, color="red")
print(g)

#Figure3c
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS")
Genecluster2<-read.table(file = "AagrOXFv2.Gene.cluster.k5.mod.all.bed",sep = "\t",header =  F)
Genecluster2$cluster<-as.factor(Genecluster2$V8)
library(dplyr)
GCLs<-as.data.frame(summary(Genecluster2$cluster))
colnames(GCLs)<-'value'
GCLs$group<-row.names(GCLs)
GCLs$ratio<-(GCLs$value/sum(GCLs$value)*100)
GCLs <- GCLs %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(GCLs$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

library(ggplot2)
ggplot(GCLs, aes(x="", y=prop, fill=group)) +theme_void()+
  geom_bar(stat="identity", width=1,color ="white") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966","gray"), labels=c("P1","P2","P3","P4","P5"))+
  geom_text(aes(y = ypos, label = round(prop,1)), color = "white", size=4) 

#Figure3d plot expression per cluster
ggplot(na.omit(merged[grep('cluster_5',merged$cluster,invert = T),]), aes(x=cluster,color=cluster,y=asinh(Ga4w_ave)))+geom_violin()+theme_classic()+
  ylab("asinh(TPM)")+
  xlab("PCG cluster")+
  scale_color_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14,vjust = 0.5), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.position = "none")+
  stat_summary(fun=median, geom="point", size=2, color="red")+
  scale_x_discrete(labels=c("P1","P2","P3","P4","P5","P6","P7","P8"))  
#Figure3e plot gene length per cluster
ggplot(na.omit(merged[grep('cluster_5',merged$cluster,invert = T),]), aes(x=cluster,color=cluster,y=log10(length)))+geom_violin()+theme_classic()+
  ylab("log10 (length)")+
  xlab("PCG cluster")+
  scale_color_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))+
  ylim(1,5)+
  theme(axis.text.y=element_text(angle=90,size=14, hjust = 0.5),
        axis.text.x=element_text(size=14,angle=90,vjust = 0.5), 
        axis.title.x=element_text(size=14,face="bold", angle=180,vjust = 2),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.position = "none")+
  stat_summary(fun=median, geom="point", size=2, color="red")+
  scale_x_discrete(labels=c("P1","P2","P3","P4"))  

#Figure3g
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS")
GoT<-read.table("GenesOverlappedbyTE25.txt", sep = "\t", header = F)
ggplot(na.omit(GoT), aes(x=V8,fill=V15))+geom_bar(color = "black")+
  theme_classic()+
  ylab("Counts")+
  xlab("PCG cluster")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33","#666666"),labels=c("T1","T2","T3","T4","T5","T6","T7","T8"))+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.title = element_text(size = 14, vjust = -0.6),
        legend.text = element_text(size = 14))+
  scale_x_discrete(labels=c("P1","P2","P3","P4","P5"))

#Figure4b
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/scale_regions")
TEcluster<- read.table("AagrOXFv2.TE.cluster.k8.mod.all.bed", header = F, sep = "\t")
TEcluster$cluster<-as.factor(TEcluster$V7)
library(dplyr)
TCLs<-as.data.frame(summary(TEcluster$cluster))
colnames(TCLs)<-'value'
TCLs$group<-row.names(TCLs)
TCLs$ratio<-(TCLs$value/sum(TCLs$value)*100)
TCLs <- TCLs %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(TCLs$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

library(ggplot2)
ggplot(TCLs, aes(x="", y=prop, fill=group)) +theme_void()+
  geom_bar(stat="identity", width=1,color ="white") +
  coord_polar("y", start=0)+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33","#666666"),labels=c("T1","T2","T3","T4","T5","T6","T7","T8"))+
  geom_text(aes(y = ypos, label = round(prop,1)), color = "white", size=4) 

#Figure4c
TEcluster<-read.table("AagrOXFv2.TE.cluster.k8.mod.bed", sep = "\t", header = F)
TEcluster$V7<-as.factor(TEcluster$V7)
TEcluster$V13<-as.factor(TEcluster$V13)
ggplot(TEcluster, aes(x=V7,color=V7,y=log10(V8)))+geom_violin()+theme_classic()+
  ylab("log10 (length)")+
  xlab("TE cluster")+
  scale_color_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33","#666666"), labels=c("T1","T2","T3","T4","T5","T6","T7","T8"))+
  ylim(1,5)+
  theme(axis.text.y=element_text(angle=90,size=14, hjust = 0.5),
        axis.text.x=element_text(size=14,angle=90,vjust = 0.5), 
        axis.title.x=element_text(size=14,face="bold", angle=180,vjust = 2),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.position = "none")+
  stat_summary(fun=median, geom="point", size=2, color="red")+
  scale_x_discrete(labels=c("T1","T2","T3","T4","T5","T6","T7","T8")) 

#Figure4e
summary(TEcluster$V7)
TEclassAll<-summary(TEcluster$V13)
TEclass1<-as.vector(summary(TEcluster[grep('cluster_1',TEcluster$V7),]$V13))
TEclass2<-as.vector(summary(TEcluster[grep('cluster_2',TEcluster$V7),]$V13))
TEclass3<-as.vector(summary(TEcluster[grep('cluster_3',TEcluster$V7),]$V13))
TEclass4<-as.vector(summary(TEcluster[grep('cluster_4',TEcluster$V7),]$V13))
TEclass5<-as.vector(summary(TEcluster[grep('cluster_5',TEcluster$V7),]$V13))
TEclass6<-as.vector(summary(TEcluster[grep('cluster_6',TEcluster$V7),]$V13))
TEclass7<-as.vector(summary(TEcluster[grep('cluster_7',TEcluster$V7),]$V13))
TEclass8<-as.vector(summary(TEcluster[grep('cluster_8',TEcluster$V7),]$V13))

TEclass<-as.data.frame(cbind(TEclassAll,TEclass1,TEclass2,TEclass3,TEclass4,TEclass5,TEclass6,TEclass7,TEclass8))
colnames(TEclass)<-c("All TEs", "cluster T1","cluster T2","cluster T3","cluster T4","cluster T5","cluster T6","cluster T7","cluster T8")
TEclass$TE_class<-row.names(TEclass)
library(reshape2)
library(ggplot2)
TEclass2<-melt(TEclass )
colnames(TEclass2)<-c("TE_class","TE_cluster", "counts")
TEclass2$counts<-as.numeric(as.character(TEclass2$counts))

ggplot(TEclass2, aes(x=TE_cluster, fill=TE_class,y=counts))+geom_bar(position="fill", stat="identity", color = "black")+theme_classic()+
  scale_fill_manual(values = c("#FF0000","#CC3366","#CC3399","#990066","#660033","#FF9900",
                               "#3399FF","#0066CC","#003366","#006699",
                               "#FFCC66","#CC9933","#996600","#CC9966","#FFCC99",
                               "#339900","#336600","#666666"))+
  ylab("ratio")+
  xlab("TE cluster")+
  labs(fill='TE class')+ 
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14,vjust = 0.5), 
        axis.title.x=element_text(size=14,face="bold",vjust = -1),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.title = element_text(size = 14, vjust = -0.6),
        legend.text = element_text(size = 14))+
  scale_x_discrete(labels=c( "All","T1","T2","T3","T4","T5","T6","T7","T8"))

#Figure4f
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS")
ToG<-read.table("TEsOverlappedbyGenes25.txt", sep = "\t", header = F)

ggplot(na.omit(ToG), aes(x=V7,fill=V21))+geom_bar(color = "black")+
  theme_classic()+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966","gray"), labels=c("P1","P2","P3","P4","P5"))+xlab("TE cluster")+
  ylab("Counts")+
  labs(fill='PCG cluster')+ 
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.title = element_text(size = 14, vjust = -0.6),
        legend.text = element_text(size = 14))+
  scale_x_discrete(labels=c("T1","T2","T3","T4","T5","T6","T7","T8"))

#Figure5
setwd("/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/reference_point_TSS")
library(ggplot2)
closeU<-read.table("ClosestUpGenes.txt",header = F,sep="\t")
closeU$V8<-as.factor(closeU$V8)
closeU$V17<-as.factor(closeU$V17)
closeU$V9<-as.factor(closeU$V9)
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
closeD$V9<-as.factor(closeD$V9)
closeD<-closeD[grep('cluster',closeD$V17),]
closeD<-closeD[grep('cluster',closeD$V8),]
closeDGenes<-closeD[grep('1',closeD$V9),]
closeDGenes1<-closeDGenes[grep('cluster_5',closeDGenes$V8,invert = T),]
closeDGenes1<-closeDGenes1[grep('cluster_5',closeDGenes1$V17,invert = T),]
closeDTEs<-closeD[grep('2',closeD$V9),]
closeDTEs1<-closeDTEs[grep('cluster_5',closeDTEs$V8,invert = T),]
closeDTEs1<-closeDTEs1[grep('cluster_8',closeDTEs1$V17,invert = T),]

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
closeTU$V9<-as.factor(closeTU$V9)
closeTD$V9<-as.factor(closeTD$V9)

UpPCG<-as.vector(summary(closeU$V9))
DownPCG<-as.vector(summary(closeD$V9))
UpTE<-as.vector(summary(closeTU$V9))
DownTE<-as.vector(summary(closeTD$V9))
ALL<-c(17223,69475)

df<-as.data.frame(cbind(ALL,UpPCG,DownPCG,UpTE,DownTE))
row.names(df)<-c('Genes','TEs')
df$category<-row.names(df)
fisher.test(df[c(1,2),c(1,2)])
fisher.test(df[c(1,2),c(1,3)])
fisher.test(df[c(1,2),c(1,4)])
fisher.test(df[c(1,2),c(1,5)])
#Figure5a
library(reshape2)
dfm<-melt(df)
ggplot(dfm,aes(x=variable,y=value,fill=category))+geom_bar(position="fill", stat="identity", color = "black")+
  scale_fill_manual(values = c("darkgreen","tomato2"), labels=c("PCG","TE"))+theme_classic()+
  ylab("Ratio")+xlab("")+labs(fill='Closest feature')+
  scale_x_discrete(labels=c("All","PCG up","PCG down","TE up","TE down"))+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14, vjust = 1, hjust = 1,angle = 45), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

#Figure5b plot closest features TE or Genes? per cluster
GCL1<-as.vector(summary(closeU[grep('cluster_1', closeU$V8),]$V9))
GCL2<-as.vector(summary(closeU[grep('cluster_2', closeU$V8),]$V9))
GCL3<-as.vector(summary(closeU[grep('cluster_3', closeU$V8),]$V9))
GCL4<-as.vector(summary(closeU[grep('cluster_4', closeU$V8),]$V9))
GCL1D<-as.vector(summary(closeD[grep('cluster_1', closeD$V8),]$V9))
GCL2D<-as.vector(summary(closeD[grep('cluster_2', closeD$V8),]$V9))
GCL3D<-as.vector(summary(closeD[grep('cluster_3', closeD$V8),]$V9))
GCL4D<-as.vector(summary(closeD[grep('cluster_4', closeD$V8),]$V9))
dfG<-as.data.frame(cbind(UpPCG,DownPCG,GCL1,GCL1D,GCL2,GCL2D,GCL3,GCL3D,GCL4,GCL4D))
row.names(dfG)<-c('PCGs','TEs')
dfG$category<-row.names(dfG)
library(reshape2)
dfGm<-melt(dfG)
ggplot(dfGm,aes(x=variable,y=value,fill=category))+geom_bar(position="fill", stat="identity", color = "black")+
  scale_fill_manual(values = c("darkgreen","tomato2"), labels=c("PCG","TE"))+theme_classic()+
  ylab("Ratio")+xlab("PCG cluster")+labs(fill='Closest feature')+
  scale_x_discrete(labels=c("All up","All down","P1 up","P1 down","P2 up","P2 down","P3 up","P3 down", "P4 up",  "P4 down"))+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14, vjust = 1, hjust = 1,angle = 45), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

#Figure5c plot closest features TE or Genes?
TCL1<-as.vector(summary(closeTU[grep('cluster_1', closeTU$V8),]$V9))
TCL2<-as.vector(summary(closeTU[grep('cluster_2', closeTU$V8),]$V9))
TCL3<-as.vector(summary(closeTU[grep('cluster_3', closeTU$V8),]$V9))
TCL4<-as.vector(summary(closeTU[grep('cluster_4', closeTU$V8),]$V9))
TCL5<-as.vector(summary(closeTU[grep('cluster_5', closeTU$V8),]$V9))
TCL6<-as.vector(summary(closeTU[grep('cluster_6', closeTU$V8),]$V9))
TCL7<-as.vector(summary(closeTU[grep('cluster_7', closeTU$V8),]$V9))
TCL1D<-as.vector(summary(closeTD[grep('cluster_1', closeTD$V8),]$V9))
TCL2D<-as.vector(summary(closeTD[grep('cluster_2', closeTD$V8),]$V9))
TCL3D<-as.vector(summary(closeTD[grep('cluster_3', closeTD$V8),]$V9))
TCL4D<-as.vector(summary(closeTD[grep('cluster_4', closeTD$V8),]$V9))
TCL5D<-as.vector(summary(closeTD[grep('cluster_5', closeTD$V8),]$V9))
TCL6D<-as.vector(summary(closeTD[grep('cluster_6', closeTD$V8),]$V9))
TCL7D<-as.vector(summary(closeTD[grep('cluster_7', closeTD$V8),]$V9))
dfT<-as.data.frame(cbind(UpTE,DownTE,TCL1,TCL1D,TCL2,TCL2D,TCL3,TCL3D,TCL4,TCL4D,TCL5,TCL5D,TCL6,TCL6D,TCL7,TCL7D))
dfT$category<-c('PCGs','TEs')
library(reshape2)
dfTm<-melt(dfT)
ggplot(dfTm,aes(x=variable,y=value,fill=category))+geom_bar(position="fill", stat="identity", color = "black")+
  scale_fill_manual(values = c("darkgreen","tomato2"), labels=c("PCG","TE"))+theme_classic()+
  ylab("Ratio")+xlab("TE cluster")+labs(fill='Closest feature')+
  scale_x_discrete(labels=c("All up","All down","T1 up","T1 down","T2 up","T2 down","T3 up", "T3 down", "T4 up", "T4 down",
                            "T5 up", "T5 down", "T6 up", "T6 down", "T7 up", "T7 down"))+
  theme(axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14, vjust = 1, hjust = 1,angle = 45), 
        axis.title.x=element_text(size=14,face="bold",hjust = 0.5),
        axis.title.y=element_text(size=14,face="bold",vjust = 2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))


#Figure5d histogram nearest features of P3
ggplot(closeUGenes1[grep('cluster_3',closeUGenes1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0), ylim = c(0,500))+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeUTEs1[grep('cluster_3',closeUTEs1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0), ylim = c(0,500))+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))


ggplot(closeDGenes1[grep('cluster_3',closeDGenes1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0,5000), ylim = c(0,500))+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeDTEs1[grep('cluster_3',closeDTEs1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0,5000), ylim = c(0,500))+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))


#Figure5 e and f histogram T4 and T5
ggplot(closeTUGenes1[grep('cluster_4',closeTUGenes1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0),ylim = c(0, 1500))+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeTUTEs1[grep('cluster_4',closeTUTEs1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0),ylim = c(0, 1500))+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))

ggplot(closeTDGenes1[grep('cluster_4',closeTDGenes1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0, 5000),ylim = c(0, 1500))+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeTDTEs1[grep('cluster_4',closeTDTEs1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0, 5000),ylim = c(0, 1500))+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))

ggplot(closeTUGenes1[grep('cluster_5',closeTUGenes1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0),ylim = c(0, 1500))+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeTUTEs1[grep('cluster_5',closeTUTEs1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(-5000,0),ylim = c(0, 1500))+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))

ggplot(closeTDGenes1[grep('cluster_5',closeTDGenes1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0, 5000),ylim = c(0, 1500))+
  xlab("Distance (bp)")+
  labs(fill='PCG cluster')+
  scale_fill_manual(values = c("#003399","#CC3333","#003300","#009966"), labels=c("P1","P2","P3","P4"))

ggplot(closeTDTEs1[grep('cluster_5',closeTDTEs1$V8),], aes(fill=V17, x=V18))+
  geom_histogram(position = "stack",binwidth = 200)+
  theme_classic()+coord_cartesian(xlim = c(0, 5000),ylim = c(0, 1500))+
  xlab("Distance (bp)")+
  labs(fill='TE cluster')+
  scale_fill_manual(values = c("#003399","#660066","#CC3333","#FF3333","#CC6666","#336633","#00CC33"),labels=c("T1","T2","T3","T4","T5","T6","T7"))
