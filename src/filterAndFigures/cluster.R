library(vcfR)
library(adegenet)
library(ggdendro)
args = commandArgs(trailingOnly=TRUE)
dataPop<-read.table(args[1],h=T,sep="\t")
gdsFile<-args[2]
popmap<-read.table(args[3],sep="\t")

colnames(popmap)<-c("sample.id","pop")

vcf<-read.vcfR(gdsFile)
genlight<-vcfR2genlight(vcf) 
nSNPs<-ncol(genlight)
             
distance<-dist(genlight,method = "euclidian")
clustering<-hclust(distance,method ="ward.D")

###########################################



dend_data<-dendro_data(clustering)
treeSegments<-dend_data$segments
treeLabels<-dend_data$labels
treeLabelsFinal<-merge(treeLabels,popmap,by.x="label",by.y="sample.id",sort=F)

splittedPath<-strsplit(gdsFile,split="/")[[1]]
max_missingString<-splittedPath[grep("max_missing",splittedPath)]
mafString<-splittedPath[grep("maf~",splittedPath)]

treeLabelsFinal$maf<-sub("^.*~","",mafString)
treeLabelsFinal$max_missing<-sub("^.*~","",max_missingString)
treeSegments$maf<-sub("^.*~","",mafString)
treeSegments$max_missing<-sub("^.*~","",max_missingString)
treeLabelsFinal$nSNPs<-nSNPs
treeLabelsFinal<-merge(treeLabelsFinal,dataPop,by="pop")
write.table(treeLabelsFinal,args[4],row.names=F,sep="\t")
write.table(treeSegments,args[5],row.names=F,sep="\t")
