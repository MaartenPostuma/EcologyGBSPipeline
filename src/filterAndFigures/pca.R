library(adegenet) 
library(vcfR)
args = commandArgs(trailingOnly=TRUE) 
dataPop<-read.table(args[1],h=T,sep="\t") 
gdsFile<-args[2] 
popmap<-read.table(args[3],sep="\t") 
vcf<-read.vcfR(gdsFile) 
genlight<-vcfR2genlight(vcf) 
colnames(popmap)<-c("sample.id","pop")
nSNPs<-ncol(genlight)
pca<-glPca(genlight,nf = nrow(genlight)-1,n.cores=4,parallel=T)
pcavar<-round(pca$eig/sum(pca$eig)*100,digits = 1) 
pcaPlot <- 
data.frame(sample.id=	row.names(pca$scores),
			PCA1 = as.data.frame(pca$scores)$PC1, # the first eigenvector 
			PCA2 = as.data.frame(pca$scores)$PC2, 
			PCA3 = as.data.frame(pca$scores)$PC3, 
			PCA1Var=pcavar[1], 
			PCA2Var=pcavar[2], 
			PCA3Var=pcavar[3], 
			nSNPs=nSNPs,# the second eigenvector 
			stringsAsFactors = FALSE) 
pcaPlot<-merge(pcaPlot,popmap,by="sample.id") 
pcavar<-pca$varpro/sum(pca$varpro,na.rm=T) 
pcaPlotFinal<-merge(pcaPlot,dataPop,by="pop",all=F) 
splittedPath<-strsplit(gdsFile,split="/")[[1]] 
max_missingString<-splittedPath[grep("max_missing",splittedPath)] 
pcaPlotFinal$max_missing<-sub("^.*~","",max_missingString) 
mafString<-splittedPath[grep("maf~",splittedPath)] 
pcaPlotFinal$maf<-sub("^.*~","",mafString) 

write.table(pcaPlotFinal,args[4],row.names=F,sep="\t")