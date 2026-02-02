args = commandArgs(trailingOnly=TRUE)
#args= c("../results/filters/fstStatsAll.tsv",
#        "../rawData/barcodes_run2.txt", 
#        "../results/filters/max_missing~1e-05/maf~0.05/populations.fst_summary.tsv",
#        "../results/filters/max_missing~1e-05/maf~0.01/populations.fst_summary.tsv")
fileList<-args[-c(1,2)]
output<-args[1] #First arg = output
barcodes<-read.table(args[2],h=T,sep="\t") #3rd arg is barcode file
listOfFiles<-list()

for(i in 1:length(fileList)){
file<-fileList[i]

fstDataTemp<-read.table(file,h=T,fill=T,sep="\t")

splittedPath<-strsplit(file,split="/")[[1]]
max_missingString<-splittedPath[grep("max_missing",splittedPath)]
mafString<-splittedPath[grep("maf~",splittedPath)]

row.names(fstDataTemp)<-fstDataTemp[,1]
fstDataTemp<-fstDataTemp[,-1]
fstDataTemp[ncol(fstDataTemp),]<-NA
row.names(fstDataTemp)[ncol(fstDataTemp)]<-colnames(fstDataTemp)[ncol(fstDataTemp)]


fstDataTemp[lower.tri(fstDataTemp)]<-t(fstDataTemp)[lower.tri(fstDataTemp)]
fstDataTemp$pop<-row.names(fstDataTemp)
fstDataTemp$maf<-sub("^.*~","",mafString)
fstDataTemp$max_missing<-sub("^.*~","",max_missingString)
listOfFiles[[i]]<-fstDataTemp
}

popStatsAll<-do.call(rbind,listOfFiles)
if(sum(colnames(barcodes)=="lat")&sum(colnames(barcodes)=="lon")){
popStatsAll<-merge(popStatsAll,unique(barcodes[,c("lat","lon","pop")]),by="pop", sort=F)

}
write.table(popStatsAll,output,row.names=F,sep="\t")
