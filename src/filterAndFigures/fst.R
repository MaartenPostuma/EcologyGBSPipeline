args = commandArgs(trailingOnly=TRUE)
fileList<-args[-c(1,2)]
output<-args[1] #First arg = output
dataBarcode<-read.table(args[2],h=T) #3rd arg is barcode file
listOfFiles<-list()

for(i in 1:length(fileList)){
file<-fileList[i]

fstDataTemp<-read.table(file,h=T,fill=T)

splittedPath<-strsplit(file,split="/")[[1]]
max_missingString<-splittedPath[grep("max_missing",splittedPath)]
mafString<-splittedPath[grep("maf~",splittedPath)]

row.names(fstDataTemp)<-fstDataTemp[,1]
fstDataTemp[,1]<-NA
fstDataTemp[ncol(fstDataTemp),]<-NA
row.names(fstDataTemp)[ncol(fstDataTemp)]<-colnames(fstDataTemp)[ncol(fstDataTemp)]

for(j in 1:nrow(fstDataTemp)){
  fstDataTemp[j,1:ncol(fstDataTemp)]<-c(rep(NA,times=j-1),fstDataTemp[j,1:(ncol(fstDataTemp)-(j-1))])
}


fstDataTemp[lower.tri(fstDataTemp)]<-t(fstDataTemp)[lower.tri(fstDataTemp)]
fstDataTemp$maf<-sub("^.*~","",mafString)
fstDataTemp$max_missing<-sub("^.*~","",max_missingString)
listOfFiles[[i]]<-fstDataTemp
}

popStatsAll<-do.call(rbind,listOfFiles)
if(sum(colnames(barcodes)=="lat")&sum(colnames(barcodes)=="lon")){
popStatsAll<-merge(unique(barcodes[,c("lat","lon","pop")]),popStatsAll,by="pop")

}
write.table(popStatsAll,output,row.names=F)