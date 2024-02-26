library(SNPRelate)

args = commandArgs(trailingOnly=TRUE)
dataPop<-read.table(args[1],h=T)
gdsFile<-args[2]
popmap<-read.table(args[3])
genofile<-snpgdsOpen(gdsFile)

colnames(popmap)<-c("sample.id","pop")

snp.id<-read.gdsn(index.gdsn(genofile,"snp.id"))
nSNPs<-length(snp.id)
pca <- snpgdsPCA(genofile, snp.id=snp.id)
pcavar<-pca$varpro/sum(pca$varpro,na.rm=T)
pcaPlot <- data.frame(sample.id = pca$sample.id,
                  PCA1 = pca$eigenvect[,1],    # the first eigenvector
                  PCA2 = pca$eigenvect[,2],
                  PCA3 = pca$eigenvect[,3],
                  PCA1Var=pcavar[1],
                  PCA2Var=pcavar[2],
                  PCA3Var=pcavar[3],
                  nSNPs=nSNps,# the second eigenvector
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