args = commandArgs(trailingOnly=TRUE)
data<-read.table(paste0(args[1]),header=F)

table2<-as.data.frame(table(data$V2))
for(pop in table2$Var1[table2$Freq>=as.numeric(args[3])]){
dataSub<-data[data$V2==pop,]
write.table(dataSub,paste0(args[2],"/",pop,"/",pop,"map.tsv"),row.names=F,col.names=F,quote=F,sep="\t")
}