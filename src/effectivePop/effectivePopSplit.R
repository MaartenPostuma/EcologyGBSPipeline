args = commandArgs(trailingOnly=TRUE)
data<-read.table(paste0(args[1]),header=F)

for(pop in unique(data$V2)){
dataSub<-data[data$V2==pop,]
write.table(dataSub,paste0(args[2],"/",pop,"/",pop,"map.tsv"),row.names=F,col.names=F,quote=F,sep="\t")
}