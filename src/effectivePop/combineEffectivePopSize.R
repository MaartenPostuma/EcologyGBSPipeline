args = commandArgs(trailingOnly=TRUE)
outFile<-args[1]
files<-args[-1]

allOut<-data.frame()
for(file in files){
test<-read.table(file,fill=T,skip = 16)



readline<-function(data,col,line){  return(  c(as.numeric(data[line,seq(col,col+3)])))}

out<-as.data.frame(t(rbind(readline(test,6,1),readline(test,4,2),
                         readline(test,4,3),readline(test,5,4),readline(test,4,5),
                         readline(test,3,7),readline(test,1,8),
                         readline(test,5,9),readline(test,1,10))))
colnames(out)<-c("meanSamples","nComparisons","Rsquared","expRsquared","estNe","paramLowCI","paramUpCI","jackKnifeLowCI","jackKnifeUpCI")    

out$maf<-c(0.05,0.02,0.01,0)
out$pop<-gsub("_Ne.txt","",gsub("^.*/","",file))
allOut<-rbind(allOut,out)
}
allOut[is.na(allOut)]<-Inf
write.table(allOut,outFile,row.names=F)