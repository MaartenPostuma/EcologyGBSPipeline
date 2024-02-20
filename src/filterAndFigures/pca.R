rm(list=ls())
library(SNPRelate)

args = commandArgs(trailingOnly=TRUE)
print(args[1])

print("test")
dataPop<-read.table(args[1],h=T)
print("test2")
gdsFile<-args[2]
popmap<-read.table(args[3])
genofile<-snpgdsOpen(gdsFile)

#dataPop<-read.table("results/stacksFiles/SNPFilterPopMap.tsv",h=T)
#popmap<-read.table("results/stacksFiles/popmapFiltered.tsv")
colnames(popmap)<-c("sample.id","pop")
#gdsFile<-"results/filters/max_missing~0.8/maf~0.0/populations.snps.gds"
genofile<-snpgdsOpen(gdsFile)

snp.id<-read.gdsn(index.gdsn(genofile,"snp.id"))
nSNPs<-length(read.gdsn(index.gdsn(genofile,"snp.id")))
pca <- snpgdsPCA(genofile, snp.id=snp.id)
pcavar<-pca$varpro/sum(pca$varpro,na.rm=T)
pcaPlot <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV1Var=pcavar[1],
                  EV2Var=pcavar[2],
                  EV3Var=pcavar[3],# the second eigenvector
                  stringsAsFactors = FALSE)
pcaPlot<-merge(pcaPlot,popmap,by="sample.id")

pcavar<-pca$varpro/sum(pca$varpro,na.rm=T)

pcaPlotFinal<-merge(pcaPlot,dataPop,by="pop",all=F)

splittedPath<-strsplit(gdsFile,split="/")[[1]]
max_missingString<-splittedPath[grep("max_missing",splittedPath)]
pcaPlotFinal$max_missing<-sub("^.*~","",max_missingString)
mafString<-splittedPath[grep("maf~",splittedPath)]
pcaPlotFinal$maf<-sub("^.*~","",mafString)
write.table(pcaPlotFinal,args[4],row.names=F,quote=F)

snpgdsClose(genofile)