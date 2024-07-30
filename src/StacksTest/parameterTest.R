library(ggplot2)
library(reshape2)
args = commandArgs(trailingOnly=TRUE)

data<-read.table(paste0(args[1],"/parameter.tsv"))
colnames(data)<-c("nVariantLoci","nLoci","coverage","M")

meltedData<-melt(data,id.vars="M")
ggplot(meltedData,aes(y=value,x=M,col=variable))+geom_line()+facet_wrap(.~variable,scales="free_y")
  ggsave(paste0(args[1],"/parameter.png"))

data<-data[order(data$M),]
diffData<-data.frame(diffVarLoci=diff(data$nLoci),diffLoci=diff(data$nVariantLoci),diffCov=diff(data$coverage),M=data$M[-1])
meltDiff<-melt(diffData,id.vars="M")
ggplot(meltDiff,aes(y=value,x=M,col=variable,group=variable))+geom_line()+facet_wrap(.~variable,scales="free_y")
  ggsave(paste0(args[1],"/parameterDiff.png"))
