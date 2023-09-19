BiocManager::install("profileplyr")
library(profileplyr)
setwd("/Volumes/groups/berger/lab/tetsuya.hisanaga/anthoceros_chromatin/data/ChIPseq/merged_bam/bigwig/scale_regions")
profileTE<-import_deepToolsMat("scale_regions1000.TEs.gz")
str(profileTE)
str(assays(profileTE))
data = cbind(assays(profileTE)[[1]],assays(profileTE)[[2]],assays(profileTE)[[3]],assays(profileTE)[[4]],assays(profileTE)[[5]],assays(profileTE)[[6]],assays(profileTE)[[7]],assays(profileTE)[[8]])
wss <- (nrow(data)-1)*sum(apply(data,2,var))

for (i in 2:12) {
  wss[i] <- sum(kmeans(data,centers=i)$withinss)
}

plot(1:11, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

profileGenes<-import_deepToolsMat("scale_regions2000.genes.gz")
str(profileGenes)
str(assays(profileGenes))
dataG = cbind(assays(profileGenes)[[1]],assays(profileGenes)[[2]],assays(profileGenes)[[3]],assays(profileGenes)[[4]],assays(profileGenes)[[5]],assays(profileGenes)[[6]],assays(profileGenes)[[7]],assays(profileGenes)[[8]])
wssG <- (nrow(dataG)-1)*sum(apply(dataG,2,var))

for (i in 2:12) {
  wssG[i] <- sum(kmeans(dataG,centers=i)$withinss)
}

plot(1:12, wssG, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

