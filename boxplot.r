library(ggplot2)
library(phytools)
args = commandArgs(trailingOnly=TRUE)
print (args)


FeaID=args[1]#paste(Dir,'3454_Lkb1_All_feature_input',sep='')
TreeID=args[2]#paste(Dir,'3454_Lkb1_All_timetree_relTimes_1_prune',sep='')
HitTaIn=args[3]

Fea=paste(FeaID,'.txt',sep='')
FeaTa=read.delim(Fea,sep=',',header=FALSE)
FeaCellOr=read.delim(paste(FeaID,'_cell.txt',sep=''),header=FALSE)
FeaGeneOr=read.delim(paste(FeaID,'_gene.txt',sep=''),header=FALSE)
row.names(FeaTa)=FeaCellOr[,1]
colnames(FeaTa)=FeaGeneOr[,1]

CelAnno=read.delim(paste(TreeID,'.txt',sep=''),sep='\t')
row.names(CelAnno)=CelAnno$Cell
FeaTaNode=merge(CelAnno,FeaTa,by='row.names')

Tree=ape::read.tree(paste(TreeID,'.nwk',sep=''))


HitTa=read.delim(HitTaIn,sep='\t')
HitGeneLs=HitTa$Gene
HitNodeLs=HitTa$Node
BooLs=HitTa$BooSco

c=1
Len=length(HitGeneLs)
print (Len)
while (c<=Len){
Gene=HitGeneLs[c]
Node=HitNodeLs[c]
Boo=BooLs[c]
if (Boo>=90){
Exp=FeaTaNode[,match(Gene, names(FeaTaNode))]
names(Exp)=FeaTaNode$Node

postscript(paste(FeaID,'_',Gene,'.ps',sep=''))
par(mfrow=c(1,2))
plotTree(Tree,mar=c(5.1,1.1,2.1,0.1))
par(mar=c(5.1,0.1,2.1,1.1))
boxplot(Exp~factor(names(Exp),levels=Tree$tip.label),horizontal=TRUE,axes=FALSE,xlim=c(1,Ntip(Tree)))
axis(1)
title(main=paste(Gene))
dev.off()
}
c=c+1
}

