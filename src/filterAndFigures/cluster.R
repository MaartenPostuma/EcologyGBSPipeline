library(SNPRelate)
library(ggdendro)
args = commandArgs(trailingOnly=TRUE)
#dataPop<-read.table("resultsBosbingelkruid/stacksFiles/SNPFilterPopMap.tsv",h=T)
dataPop<-read.table(args[1],h=T)
#gdsFile<-"resultsBosbingelkruid/filters/max_missing~0.8/maf~0.025/populations.snps.gds"
gdsFile<-args[2]
popmap<-read.table(args[3])
#popmap<-read.table("resultsBosbingelkruid/stacksFiles/popmapFiltered.tsv")
genofile<-snpgdsOpen(gdsFile)

colnames(popmap)<-c("sample.id","pop")

snp.id<-read.gdsn(index.gdsn(genofile,"snp.id"))
nSNPs<-length(snp.id)
diss<-snpgdsDiss(genofile,autosome.only=F)

#This might have to do to with the practice dataSet will remove for real thing
#checkForbrokenInds<-which(is.na(diss$diss))[1]
#diss$diss<-diss$diss[-checkForbrokenInds,-checkForbrokenInds]
#diss$sample.id<-diss$sample[-checkForbrokenInds]
#checkForBrokenInds2<-which(!complete.cases(diss$diss))
#diss$diss<-diss$diss[-checkForBrokenInds2,-checkForBrokenInds2]
#diss$sample.id<-diss$sample[-checkForBrokenInds2]
########################################
splittedPath<-strsplit(gdsFile,split="/")[[1]]
max_missingString<-splittedPath[grep("max_missing",splittedPath)]
mafString<-splittedPath[grep("maf~",splittedPath)]


tryCatch({
hc<-snpgdsHCluster(diss)

fit <- hclust(as.dist(hc$dist), method="ward.D") 
dend_data <- dendro_data(fit, type = "rectangle")
treeSegments<-dend_data$segments
treeLabels<-dend_data$labels
treeLabelsFinal<-merge(treeLabels,popmap,by.x="label",by.y="sample.id",sort=F)


treeLabelsFinal$maf<-sub("^.*~","",mafString)
treeLabelsFinal$max_missing<-sub("^.*~","",max_missingString)
treeSegments$maf<-sub("^.*~","",mafString)
treeSegments$max_missing<-sub("^.*~","",max_missingString)
treeLabelsFinal$nSNPs<-nSNPs
treeLabelsFinal<-merge(treeLabelsFinal,dataPop,by="pop")
write.table(treeLabelsFinal,args[4],row.names=F,quote=F)
write.table(treeSegments,args[5],row.names=F,quote=F)

snpgdsClose(genofile)
},error=function(cond) {
            message(conditionMessage(cond))
            write.table(data.frame(pop=NA,label=NA,x=NA,y=NA,max_missing=max_missingString,maf=mafString, nSNPs=0, metaPop=NA),args[4],row.names=F,quote=F)
            write.table(data.frame(x=NA,y=NA,xend=NA,yend=NA,max_missing=max_missingString,maf=mafString),args[5],row.names=F,quote=F)

})