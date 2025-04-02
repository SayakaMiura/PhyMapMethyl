args = commandArgs(trailingOnly=TRUE)
print (args)
HitTaIn=args[1]
HitTa=read.delim(HitTaIn,sep='\t')
NodeLs=unique(HitTa$Node)
for (Node in NodeLs){
HitTa1=HitTa[HitTa$Node==Node,]
if (length(HitTa1$Gene)==1){
HitTa1Rep=HitTa1
HitTa1Rep$Gene=c('duplicate')
HitTa1=rbind(HitTa1,HitTa1Rep)
}
row.names(HitTa1)=HitTa1$Gene
A=as.matrix(HitTa1[,5:(ncol(HitTa))])
png(paste(HitTaIn,Node,'.png',sep=''))
heatmap(A,main=paste(Node))
dev.off()
}
