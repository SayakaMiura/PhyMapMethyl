from Bio import Phylo
from Bio.Phylo.Consensus import *
from io import StringIO
import numpy as np
import os
import glob
import sys

class Tree(object):
    def __init__(self):
        self.message={}
    def ReadNwkBra(self,Tree):
      #  print (Tree)
     #   open('a','r').readlines()
        tree = Phylo.read(Tree, "newick")
        C=1
        for i in tree.find_clades():
           if i.name==None:
                i.name=C
       #print (i.name)
           C+=1
        Tips=tree.get_terminals()
      #  print (Tips)
      #  open('a','r').readlines()
        Dec2Anc={}
        Bra2Len={}
        for Tip in Tips:  
           Root2TipLs=tree.get_path(Tip.name)
          # print (Tip,Root2TipLs)
         #  open('a','r').readlines()
           D=1
           Len=len(Root2TipLs)
          # print (Tip.name,Root2TipLs)
           if Len==1:
           #    print ('h',Root2TipLs[0].name)
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
      #  print (Dec2Anc)
      #  open('a','r').readlines()             
        return Dec2Anc,Tips1,Bra2Len
    
    def InvertDic1(self,St2Seq):
        Hap2ID={}
        for St in St2Seq:
            Seq=St2Seq[St]
            Hap2ID[Seq]=Hap2ID.get(Seq,[])+[St]
        return Hap2ID    
    def GetClade(self,Anc,Dec2Anc,TipLs):
        In=[]
        if Anc in TipLs: In.append(Anc)
      #  print (Dec2Anc)
        for i in TipLs:
           # print (i,Dec2Anc[str(i).strip()])
            F='n'
            A=Dec2Anc.get(str(i).strip(),'')
           # print (A,Anc)
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
    def Che(self,FeaLs,Clade):
        Good=[]
        Oth=[]
        for i in Clade:
            if i in FeaLs: Good.append(i)
            else: Oth.append(i)
        Good.sort()
        print (len(Good),(len(FeaLs)))
        if len(Good)>=(len(FeaLs)-1): Good=[]
        #print (FeaLs)
        return Good,Oth  
        
            
        ###
    def readCellAnno(self,CellAnno):
        Node2CellLs={}
        CellAnno=open(CellAnno,'r').readlines()[1:]
        for i in CellAnno:
            i=i.split('\t')
            Cell=i[0].strip()
            Node=i[1].strip()
            Node2CellLs[Node]=Node2CellLs.get(Node,[])+[Cell]
        return Node2CellLs  
    def GetCellFromNode(self,Dic,NodeLs):
            CellLs=[]
            for Node in NodeLs:
                 CellLs+=Dic.get(Node,[])
            return CellLs                
    def MakeESLinFromTree(self,NWK,CellAnno):                
       # NWK=sys.argv[1]#Dir+ID+'.nwk'
        Node2CellLs=self.readCellAnno(CellAnno)
      #  print (len(Node2CellLs))
       # Dir=sys.argv[2]#NWK.replace(ID+'.nwk','')#'/home/sayaka/Desktop/ESL'+os.sep #change
       # ID=Dir.split(os.sep)[-1]#NWK.split(ps.sep)[-1][:-4]#'C551' #change
        Out=NWK[:-4]+'_ingroup.txt'
       # print (ID,Dir,NWK)
        
      #  FeaFaLs=glob.glob(Dir+os.sep+ID+'-*.txt')
      #  print (len(FeaFaLs))
      #  Dir=Dir[:(-1*(len(ID)))]
      #  print (Dir)
      #  #open('a','r').readlines()
        
           
        out=['\t'.join(['ID','Set1','Set2','BranchLength'])+'\n']
        Dec2Anc,Tips,Bra2Len= self.ReadNwkBra(NWK)
      #  print (Dec2Anc,len(Tips))
        #print (Bra2Len)
        #print (len(Tips))
      #  FeaLs=[]
      #  for i in FeaFaLs:
      #      SampID=i.replace(Dir+ID+os.sep,'')[:-4]
      #      if SampID in Tips:
      #          FeaLs.append(SampID)
      #      else: print ('sample not found in tree ',SampID)        
      #  print ('feature','tree',len(FeaLs),len(Tips))
        #open('a','r').readlines()
        Anc2Dec=self.InvertDic1(Dec2Anc)
     #   print (Anc2Dec)
        FeaLs=Tips
        Done=[]
        #print (Anc2Dec)
        NodeMap=['ID\tNode\tAncestor\tSet1\tSet2\n']
        for Anc in Dec2Anc:

            All,Neg=self.GetClade(Anc,Dec2Anc,FeaLs)
            #print (Anc,All,Neg)
            NodeID=''
           # All,Neg=Che(FeaLs,Clade)
          #  print (All)
            if len(All)>0 and len(Neg)>0:# and Done.count(All)==0: 
                NodeID='Node'+str(Anc)
                Done.append(All)
  
                AA=Dec2Anc.get(Anc,'root')
                NodeMap.append('\t'.join(map(str,[NodeID,Anc,AA]))+'\t'+';'.join(All)+'\t'+';'.join(Neg)+'\n')
           # print (Anc,AA,Dec2Anc)
              #  print ('get cell')
               # Allcell=self.GetCellFromNode(Node2CellLs,All)
               # Negcell=self.GetCellFromNode(Node2CellLs,Neg)
               # print (len(Allcell),len(Negcell))
                out.append(NodeID+'\t'+';'.join(All)+'\t'+';'.join(Neg)+'\n')
            
       # print ('out1')    
        OutF=open(Out,'w')
        OutF.write(''.join(out))
        OutF.close()
      #  print ('out1')
        OutF=open(Out[:-4]+'NodeMap.txt','w')
        OutF.write(''.join(NodeMap))
        OutF.close() 
       # Cell2Node=self.InvertDic1(Node2CellLs)
        return Node2CellLs       
