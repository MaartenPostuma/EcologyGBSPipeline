args = commandArgs(trailingOnly=TRUE)

fileList<-args[-c(1,2)]
dataPop<-read.table(args[2],h=T)
output<-args[1]


listOfFiles<-list()
for(i in 1:length(fileList)){
file<-fileList[i]

popDataTemp<-read.table(file,h=T)

splittedPath<-strsplit(file,split="/")[[1]]
max_missingString<-splittedPath[grep("max_missing",splittedPath)]
mafString<-splittedPath[grep("maf~",splittedPath)]

popDataTemp<-read.table(file,h=T)
popDataTemp$max_missing<-sub("^.*~","",max_missingString)
popDataTemp$maf<-sub("^.*~","",mafString)
popDataTemp$percPoly<-(popDataTemp$Polymorphic_Sites/popDataTemp$Variant_Sites)*100
listOfFiles[[i]]<-read.table(file,h=T)
}

popStatsAll<-do.call(rbind,listOfFiles)
write.table(popStatsAll,output,row.names=F)