rm(list=ls())

#TODO make the SNPFIlter popmap but first rewrite most of the SNPFiltering scripts!
args = commandArgs(trailingOnly=TRUE)


barcodes<-read.table(args[1],h=T,sep="\t")
barcodes$run<-sub("_R1.fq.gz","",barcodes$rawR1)

barcodeStacks<-data.frame(sample=barcodes$sample,
                          barcode1=paste0(barcodes$barcode1,"C"),
                          barcode2=paste0(barcodes$barcode2,"C"))
write.table(barcodeStacks,paste0(args[2],"/barcodeStacks.tsv"),row.names = F,col.names = F,quote=F,sep="\t")

for(run in barcodes$run){
  barcodesSub<-barcodes[barcodes$run==run,]
  barcodeStacks<-data.frame(barcode1=paste0(barcodesSub$barcode1,"C"),
                            barcode2=paste0(barcodesSub$barcode2,"C"),
                            sample=barcodesSub$sample)
  barcodeStacks$sample[duplicated(barcodeStacks$sample)]<-paste0(barcodeStacks$sample[duplicated(barcodeStacks$sample)],"DUP")
  write.table(barcodeStacks,paste0(args[2],"/barcodeStacks",run,".tsv"),row.names = F,col.names = F,quote=F,sep="\t")
  popmapStacksSub<-data.frame(sample=barcodes$sample[duplicated(barcodes$sample)==F],pop=barcodes$pop[duplicated(barcodes$sample)==F])
  write.table(popmapStacksSub,paste0(args[2],"/popmap",run,".tsv"),row.names = F,col.names = F,quote=F,sep="\t")    
}


popmapStacks<-data.frame(sample=barcodes$sample[duplicated(barcodes$sample)==F],pop=barcodes$pop[duplicated(barcodes$sample)==F])
write.table(popmapStacks,paste0(args[2],"/popmap.tsv"),row.names = F,col.names = F,quote=F,sep="\t")


popmapSNPFilter<-data.frame(pop=unique(barcodes$pop),metaPop=barcodes$metaPop[duplicated(barcodes$pop)==F])
write.table(popmapSNPFilter,paste0(args[2],"/SNPFilterPopMap.tsv"),row.names =F,col.names =T,quote=F,sep="\t")

