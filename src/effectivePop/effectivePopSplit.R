args = commandArgs(trailingOnly=TRUE)
data<-read.table(paste0(args[1]))

for(pop in data$V2){
dataSub<-data[data$V2==pop,]
write.table(dataSub,paste0(args[2],"/",pop,"map.tsv"),row.names=F,colnames=F,quote=F,sep="\t")
}