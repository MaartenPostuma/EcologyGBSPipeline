library(SNPRelate)

args = commandArgs(trailingOnly=TRUE)
dataPop<-read.table(args[1],h=T)
gdsFile<-args[2]
popmap<-read.table(args[3])
genofile<-snpgdsOpen(gdsFile)

colnames(popmap)<-c("sample.id","pop")

snp.id<-read.gdsn(index.gdsn(genofile,"snp.id"))
nSNPs<-length(snp.id)

diss<-snpgdsDiss(genofile)
hc<-snpgdsHCluster(diss)

fit <- hclust(as.dist(hc$dist), method="ward.D") 
dend_data <- dendro_data(fit, type = "rectangle")
treeSegments<-dend_data$segments
treeLabels<-dend_data$labels
treeLabelsFinal<-merge(treeLabels,popmap,by.x="label",by.y="V1",sort=F)


splittedPath<-strsplit(gdsFile,split="/")[[1]]
max_missingString<-splittedPath[grep("max_missing",splittedPath)]
treeLabelsFinal$max_missing<-sub("^.*~","",max_missingString)
treeSegments$max_missing<-sub("^.*~","",max_missingString)
mafString<-splittedPath[grep("maf~",splittedPath)]
treeLabelsFinal$maf<-sub("^.*~","",mafString)
treeSegmentsl$maf<-sub("^.*~","",mafString)
treeLabelsFinal$nSNPs<-nSNPs
write.table(treeLabelsFinal,args[4],row.names=F,quote=F)
write.table(treeSegments,args[5],row.names=F,quote=F)

snpgdsClose(genofile)