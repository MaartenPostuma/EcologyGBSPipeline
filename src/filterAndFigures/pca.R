rm(list=ls())
library(SNPRelate)

args = commandArgs(trailingOnly=TRUE)



dataPop<-read.table(args[1])
genofile<-snpgdsOpen(args[2])

#dataPop<-read.table("results/stacksFiles/SNPFilterPopMap.tsv")
gdsFile<-"results/filters/max_missing~0.8/maf~0.0/populations.snps.gds"
genofile<-snpgdsOpen("results/filters/max_missing~0.8/maf~0.0/populations.snps.gds")


ranges<-sub('vcf/filter1//filter1_','',sub('.recode.vcf','',fileList))

tryCatch(for(i in fileList){
range<-sub('vcf/filter1//filter1_','',sub('.recode.vcf','',i))
vcfPath<-i
dataLoc<-read.table(dataPop,h=T)#Needs to be taken from config BLEGH
if(!file.exists(sub(".vcf",".gds",vcfPath))){
invisible(snpgdsVCF2GDS(vcfPath, sub(".vcf",".gds",vcfPath), method="biallelic.only"))}
invisible(genofile <- snpgdsOpen(sub(".vcf",".gds",vcfPath)))
snp.id<-read.gdsn(index.gdsn(genofile,"snp.id"))
nSNPs<-length(read.gdsn(index.gdsn(genofile,"snp.id")))
invisible(pca <- snpgdsPCA(genofile, snp.id=snp.id))
pcaPlot <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],# the second eigenvector
                  stringsAsFactors = FALSE)
pcaPlot$pop<-as.character(sub('_.*$', '', pcaPlot$sample.id))
pcavar<-pca$varpro/sum(pca$varpro,na.rm=T)

pcaPlotFinal<-merge(pcaPlot,dataLoc,by="pop",all=F)
print(ggplot(pcaPlotFinal,aes(x=EV1,y=EV2,col=metaPop,label=sample.id))+geom_text()+
theme_light()+ylab(paste("PCA 2 ",round(pcavar[2],3)*100,"%",sep=""))+ggtitle(vcfPath)+xlab(paste("PCA 1 ",round(pcavar[1],3)*100,"%",sep=""))+scale_colour_discrete("meta populations")+
ggtitle(paste("number of SNPs =",nSNPs,"\n maxMissing =",range)))
snpgdsClose(genofile)

})
