library(l1ou)
library(ape)
library(phytools) 
args = commandArgs(trailingOnly=TRUE)
print (args)
ESL=args[1]
Tree=args[2]
print (Tree)
Tr0=ape::read.tree(Tree)
Max=1
print (Max)
Tr=force.ultrametric(Tr0, method="extend")
HitS=4 ###
#BooT=100
HitE=HitS+length(Tr$tip.label)-1
Hit=read.delim(ESL)

AddedGene=0
Gene=1
TotGene=length(Hit$Gene)
#TotGene=10
while (Gene<=TotGene){ 
Sub= Hit[Gene, HitS:HitE]

In=adjust_data(Tr, t(Sub), normalize = TRUE, quietly = TRUE)

eModel=c()
result=c()
tryCatch({

eModel <- estimate_shift_configuration(In$tree, In$Y,criterion="BIC",max.nShifts = Max,quietly = TRUE)},error=function(e){}) 


Ed=eModel$tree$edge

if (is.null(Ed)==FALSE){
Res=c(Hit[Gene,1],eModel$alpha,eModel$sigma2,eModel$shift.values,eModel$alpha*eModel$shift.values)
print (Res)
if (length(Res)==5){
if (abs(eModel$alpha*eModel$shift.values)>0.2 & eModel$sigma2<0.1){
colnames(Ed)=c('Anc','ID')
TipLa=data.frame(eModel$tree[4], seq(1,length(eModel$tree$tip.label),by=1))
colnames(TipLa)=c('TipLabel','ID')
A=merge(Ed,TipLa,by='ID',all=TRUE)
if (AddedGene==0){

Shi=t(rbind(eModel$shift.configuration,t(eModel$alpha*eModel$shift.values),rep(Hit[Gene,1],length(eModel$shift.values))))
colnames(Shi)=c('branch','S.v','Gene')
DecAnc=cbind(A,rep(Hit[Gene,1],length(A$ID)))
colnames(DecAnc)=c('ID','Anc','TipLabel','Gene')
Node=data.frame(rep(Hit[Gene,1],length(eModel$tree$tip.label)),eModel$tree$tip.label,eModel$l1ou.options$Z)
AddedGene=AddedGene+1
} else{

ShiAdd=t(rbind(eModel$shift.configuration,t(eModel$alpha*eModel$shift.values),rep(Hit[Gene,1],length(eModel$shift.values))))
colnames(Shi)=c('branch','S.v','Gene')
Shi=rbind(Shi,ShiAdd)
DecAncAdd=cbind(A,rep(Hit[Gene,1],length(A$ID)))
colnames(DecAncAdd)=c('ID','Anc','TipLabel','Gene')
DecAnc=rbind(DecAnc,DecAncAdd)
NodeAdd=data.frame(rep(Hit[Gene,1],length(eModel$tree$tip.label)),eModel$tree$tip.label,eModel$l1ou.options$Z)
Node=rbind(Node,NodeAdd)
}}}}
Gene=Gene+1
}

write.table(Shi,paste(ESL,'_Shift.txt',sep=''),quote=FALSE,sep='\t', row.names = FALSE)

write.table(Node,paste(ESL,'_TreeNode.txt',sep=''),quote=FALSE,sep='\t', row.names = FALSE)

