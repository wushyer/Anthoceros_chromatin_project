a=read.table("curated.info",head=T,sep="\t")
head(a)
library(tidyverse)
library(patchwork)
library(ggridges)
library(ggsci)


methy=ggplot(a, aes(Methylation,chr,fill=chr), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5)+scale_fill_npg() +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K27me1=ggplot(a, aes(H3K27me1,chr,fill=chr), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5)+scale_fill_npg() +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )

H3K27me3=ggplot(a, aes(H3K27me3,chr,fill=chr), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5)+scale_fill_npg() +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K36me3=ggplot(a, aes(H3K36me3,chr,fill=chr), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5)+scale_fill_npg() +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K4me3=ggplot(a, aes(H3K4me3,chr,fill=chr), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5)+scale_fill_npg() +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K9me1=ggplot(a, aes(H3K9me1,chr,fill=chr), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5)+scale_fill_npg() +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )

methy+H3K9me1+H3K27me1+H3K27me3+H3K4me3+H3K36me3

###### volin Compartment####
methy.Compartment=ggplot(a, aes(Methylation,Compartment,fill=Compartment), guide="none") +
  geom_violin(aes(fill=Compartment),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    #axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
    plot.caption = element_text(hjust = 0)
  )

H3K27me1.Compartment=ggplot(a, aes(H3K27me1,Compartment,fill=Compartment), guide="none") +
  geom_violin(aes(fill=Compartment),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    #axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
    plot.caption = element_text(hjust = 0)
  )

H3K27me3.Compartment=ggplot(a, aes(H3K27me3,Compartment,fill=Compartment), guide="none") +
  geom_violin(aes(fill=Compartment),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K36me3.Compartment=ggplot(a, aes(H3K36me3,Compartment,fill=Compartment), guide="none") +
  geom_violin(aes(fill=Compartment),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K4me3.Compartment=ggplot(a, aes(H3K4me3,Compartment,fill=Compartment), guide="none") +
  geom_violin(aes(fill=Compartment),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K9me1.Compartment=ggplot(a, aes(H3K9me1,Compartment,fill=Compartment), guide="none") +
  geom_violin(aes(fill=Compartment),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    #axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
methy.Compartment+H3K9me1.Compartment+H3K27me1.Compartment+H3K27me3.Compartment+H3K4me3.Compartment+H3K36me3.Compartment


####volin TAD#####
a=read.table("TAD.curated.info",head=T,sep="\t")
methy.TAD=ggplot(a, aes(Methylation,TAD,fill=TAD), guide="none") +
  geom_violin(aes(fill=TAD),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )

H3K27me1.TAD=ggplot(a, aes(H3K27me1,TAD,fill=TAD), guide="none") +
  geom_violin(aes(fill=TAD),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )

H3K27me3.TAD=ggplot(a, aes(H3K27me3,TAD,fill=TAD), guide="none") +
  geom_violin(aes(fill=TAD),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K36me3.TAD=ggplot(a, aes(H3K36me3,TAD,fill=TAD), guide="none") +
  geom_violin(aes(fill=TAD),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K4me3.TAD=ggplot(a, aes(H3K4me3,TAD,fill=TAD), guide="none") +
  geom_violin(aes(fill=TAD),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
H3K9me1.TAD=ggplot(a, aes(H3K9me1,TAD,fill=TAD), guide="none") +
  geom_violin(aes(fill=TAD),draw_quantiles = c(0.5))+coord_flip()+scale_fill_npg()  +
  theme_classic()+
  theme(
    text = element_text(size=12, color = "black", face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.caption = element_text(hjust = 0)
  )
methy.TAD+H3K9me1.TAD+H3K27me1.TAD+H3K27me3.TAD+H3K4me3.TAD+H3K36me3.TAD
