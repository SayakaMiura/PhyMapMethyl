from Bio import Phylo
from Bio.Phylo.Consensus import *
from io import StringIO
import numpy as np
import os
import glob
import sys


NWK=sys.argv[1]
Dir=sys.argv[2]#NWK.replace(ID+'.nwk','')
ID=sys.argv[3]#NWK.split(ps.sep)[-1][:-4]
Out=NWK[:-4]+'.txt'
print (ID,Dir,NWK)

FeaFaLs=glob.glob(Dir+os.sep+ID+'-*.txt')
print ('sample count',len(FeaFaLs))

def ReadNwkBra(Tree):
    tree = Phylo.read(Tree, "newick")
    C=1
    for i in tree.find_clades():
       if i.name==None:
            i.name=C
   #print (i.name)
       C+=1
    Tips=tree.get_terminals()
    Dec2Anc={}
    Bra2Len={}
    for Tip in Tips:  
       Root2TipLs=tree.get_path(Tip.name)
       D=1
       Len=len(Root2TipLs)

       if Len==1:

           Dec2Anc[Root2TipLs[0].name]='root'
           Bra2Len[Root2TipLs[0].name]=Root2TipLs[0].branch_length
   
       while D<Len:
           Dec2Anc[Root2TipLs[D].name]=Root2TipLs[D-1].name
           Bra2Len[Root2TipLs[D].name]=Root2TipLs[D].branch_length
           D+=1
       Dec2Anc[Root2TipLs[0].name]='root'
       Bra2Len[Root2TipLs[0].name]=Root2TipLs[0].branch_length
    Tips1=[]
    for i in Tips:
        Tips1.append(i.name)             
    return Dec2Anc,Tips1,Bra2Len

def InvertDic1(St2Seq):
    Hap2ID={}
    for St in St2Seq:
        Seq=St2Seq[St]
        Hap2ID[Seq]=Hap2ID.get(Seq,[])+[St]
    return Hap2ID    
def GetClade(Anc,Dec2Anc,TipLs):
    In=[]

    for i in TipLs:

        F='n'
        A=Dec2Anc.get(str(i).strip(),'')

        if A ==Anc: In.append(i)
        else:
            while A in Dec2Anc:
                A=Dec2Anc[A]
                if A ==Anc: 
                    In.append(i)
                    break
    Oth=[]
    for i in TipLs:
       if i not in In: Oth.append(i)                 
    return In,Oth                
def Che(FeaLs,Clade):
    Good=[]
    Oth=[]
    for i in Clade:
        if i in FeaLs: Good.append(i)
        else: Oth.append(i)
    Good.sort()
    print (len(Good),(len(FeaLs)))
    if len(Good)>=(len(FeaLs)-1): Good=[]

    return Good,Oth 
def AddID(Ls,Add):
   A=[]
   for i in Ls:
       A.append(Add+'-'+i)
   return A         
out=['\t'.join(['ID','Node','Anc','Set1','Set2','BranchLength','PhyMapNodeID'])+'\n']
out1=['\t'.join(['ID','Set1','Set2'])+'\n']
Dec2Anc,Tips,Bra2Len= ReadNwkBra(NWK)

FeaLs=[]
for i in FeaFaLs:
    SampID=i.replace(Dir+os.sep+ID+'-','')[:-4]

    if SampID in Tips:
        FeaLs.append(SampID)
    elif SampID !='N': print ('sample not found in tree ',SampID)        
print ('feature','tree',len(FeaLs),len(Tips))

Anc2Dec=InvertDic1(Dec2Anc)
Done=[]

for Anc in Dec2Anc:
    All,Neg=GetClade(Anc,Dec2Anc,FeaLs)

    NodeID=''

    if len(All)>1 and len(Neg)>1 and Done.count(All)==0: 
        NodeID='Node'+str(Anc)
        Done.append(All)
    else: 
        All=[]
        Neg=[]    
 
    AA=Dec2Anc.get(Anc,'root')

    if len(All)>1:
        All.sort()

    out.append('\t'.join(map(str,[NodeID,Anc,AA]))+'\t'+';'.join(All)+'\t'+';'.join(Neg)+'\t'+str(Bra2Len[Anc])+'\t'+''.join(All)+'\n')
    if NodeID!='':
        All1=AddID(All,ID)
        Neg1=AddID(Neg,ID)
        Neg1.append(ID+'-N')
        Neg1.append(ID+'-N1')
        out1.append(NodeID+'\t'+';'.join(All1)+'\t'+';'.join(Neg1)+'\n')
Tip1=AddID(FeaLs,ID)        
out1.append('Root\t'+';'.join(Tip1)+'\t'+';'.join([ID+'-N',ID+'-N1'])+'\n')
out1.append('Outgroup\t'+ID+'-N\n')#;'+ID+'-N1\n')        
#print (out)    
OutF=open(Out,'w')
OutF.write(''.join(out1))
OutF.close()
OutF=open(Out[:-4]+'_NodeID.txt','w')
OutF.write(''.join(out)) 
OutF.close()
