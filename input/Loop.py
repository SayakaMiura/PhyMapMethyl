import numpy as np
import os
import glob
import sys
import pandas as pd
import shutil

class Loop(object):
    def __init__(self):
        self.message={}
    def UpInput0(self,Fea,GeneIn,AllHit,BooCut,Rep):
       # SuppTa=pd.read_csv(Supp,sep='\t')
       # print (SuppTa)
       # Len=len(SuppTa)
       # c=0
       # SelGeneLs=[]
       # while c<Len:
       #    # print (list(SuppTa['Gene'])[c])
       #     Gene=list(SuppTa['Gene'])[c]+'\t'+str(list(SuppTa['Position from 0'])[c])
       #     Boo=float(list(SuppTa['BooSco'])[c])
       #    # print (Gene,Boo)
       #     if Boo>=BooCut:
       #          SelGeneLs.append(Gene)
       #     c+=1
       # if len(SelGeneLs)==0: return 'Done','Done','Done'   
       # else:
            Dir=GeneIn.replace(GeneIn.split(os.sep)[-1],'')
            if os.path.exists(Dir+Rep)!=True: os.mkdir(Dir+Rep)
            UpFea=Dir+Rep+os.sep+Rep+'_feature_input.txt'
            UpCell=Dir+Rep+os.sep+Rep+'_feature_input_cell.txt'
            UpGene=Dir+Rep+os.sep+Rep+'_feature_input_gene.txt'
            if os.path.exists(GeneIn[:-9]+'_cell.txt')==True:shutil.copy2(GeneIn[:-9]+'_cell.txt',UpCell)
           # print (Dir+Rep)
            if os.path.exists(UpFea)!=True:
                HitTa=pd.read_csv(AllHit,sep='\t')
                HitPos=list(HitTa['Position from 0'])
          #  print (Fea)
                print ('prune hit',len(HitPos),HitPos[:10])
                Fea1=Fea.drop(map(int,HitPos),axis=1)
                print (Fea1)
                Fea1.to_csv(UpFea,header=False,index=False)   
                GeneLs=open(GeneIn,'r').readlines()
                GeneLs1=np.array(GeneLs)
                GeneLs2=np.delete(GeneLs1,HitPos)
                OutF=open(UpGene,'w')
                OutF.write(''.join(list(GeneLs2)))
                OutF.close()
            return UpFea,UpCell,UpGene         
    def UpInput_copy(self,Fea,GeneIn,AllHit,BooCut,Rep):
            Dir=GeneIn.replace(GeneIn.split(os.sep)[-1],'')
            if os.path.exists(Dir+Rep)!=True: os.mkdir(Dir+Rep)
            UpFea=Dir+Rep+os.sep+Rep+'_feature_input.txt'
            UpCell=Dir+Rep+os.sep+Rep+'_feature_input_cell.txt'
            UpGene=Dir+Rep+os.sep+Rep+'_feature_input_gene.txt'
            shutil.copy2(GeneIn[:-9]+'_gene.txt',UpGene)
            if os.path.exists(GeneIn[:-9]+'_cell.txt')==True: shutil.copy2(GeneIn[:-9]+'_cell.txt',UpCell)
            if os.path.exists(UpFea)!=True:
                shutil.copy2(GeneIn[:-9]+'.txt',UpFea)
            print ('too many hit, so input was not updated')
         #   print (Dir+Rep)
         #   HitTa=pd.read_csv(AllHit,sep='\t')
         #   HitPos=list(HitTa['Position from 0'])
         #   print (Fea)
         #   print ('prune hit',len(HitPos),HitPos[:10])
         #   Fea1=Fea.drop(map(int,HitPos),axis=1)
         #   print (Fea1)
         #   Fea1.to_csv(UpFea,header=False,index=False)   
         #   GeneLs=open(GeneIn,'r').readlines()
         #   GeneLs1=np.array(GeneLs)
         #   GeneLs2=np.delete(GeneLs1,HitPos)
         #   OutF=open(UpGene,'w')
         #   OutF.write(''.join(list(GeneLs2)))
         #   OutF.close()
            return UpFea,UpCell,UpGene               
    def UpInput(self,Fea,GeneIn,AllHit,Supp,BooCut,Rep):
        SuppTa=pd.read_csv(Supp,sep='\t')
        print (SuppTa)
        Len=len(SuppTa)
        c=0
        SelGeneLs=[]
        while c<Len:
           # print (list(SuppTa['Gene'])[c])
            Gene=list(SuppTa['Gene'])[c]+'\t'+str(list(SuppTa['Position from 0'])[c])
            Boo=float(list(SuppTa['BooSco'])[c])
           # print (Gene,Boo)
            if Boo>=BooCut:
                 SelGeneLs.append(Gene)
            c+=1
        if len(SelGeneLs)==0: return 'Done','Done','Done'   
        else:
            Dir=GeneIn.replace(GeneIn.split(os.sep)[-1],'')
            if os.path.exists(Dir+Rep)!=True: os.mkdir(Dir+Rep)
            UpFea=Dir+Rep+os.sep+Rep+'_feature_input.txt'
            UpCell=Dir+Rep+os.sep+Rep+'_feature_input_cell.txt'
            UpGene=Dir+Rep+os.sep+Rep+'_feature_input_gene.txt'
            if os.path.exists(GeneIn[:-9]+'_cell.txt')==True: shutil.copy2(GeneIn[:-9]+'_cell.txt',UpCell)
            print (Dir+Rep)
            HitTa=pd.read_csv(AllHit,sep='\t')
            HitPos=list(HitTa['Position from 0'])
            print (Fea)
            print ('prune hit',len(HitPos),HitPos[:10])
            Fea1=Fea.drop(map(int,HitPos),axis=1)
            print (Fea1)
            Fea1.to_csv(UpFea,header=False,index=False)   
            GeneLs=open(GeneIn,'r').readlines()
            GeneLs1=np.array(GeneLs)
            GeneLs2=np.delete(GeneLs1,HitPos)
            OutF=open(UpGene,'w')
            OutF.write(''.join(list(GeneLs2)))
            OutF.close()
            return UpFea,UpCell,UpGene 
    def UpInputBest(self,Fea,GeneIn,AllHit,BetaLs,BooCut,Rep):
      c=0
      Dir=GeneIn.replace(GeneIn.split(os.sep)[-1],'')
      for Supp in BetaLs:
       # print (Supp)
        SuppTa=pd.read_csv(Supp,sep='\t',header=None)
        SuppTa_s = SuppTa.sort_values(2)
        Small = SuppTa_s.head(1)
        Large = SuppTa_s.tail(1)
        
      #  print (Small,Large)
        if c==0:
           BestAll= pd.concat([Small, Large])
           c+=1
        else:
           Add = pd.concat([Small, Large])      
           BestAll= pd.concat([BestAll, Add])
      print (BestAll)
      abs_column = BestAll[2].abs()
      abs_df = pd.DataFrame({'Abs': abs_column})
      combined_df0 = pd.concat([BestAll, abs_df], axis=1)
      combined_df = combined_df0.dropna()
      combined_df.to_csv(Dir+'allBeta.txt',sep='\t',index=False)
      combined_df_s = combined_df.sort_values('Abs')
      RmLine = combined_df_s.tail(1)
      print (RmLine)
      HitPos=[list(RmLine[0])[0]-1]
      print ('prune best',HitPos,RmLine)
      OutF=open(Dir+'rmFeature.txt','w')
      OutF.write(str(HitPos)+'\n'+str(RmLine))
      OutF.close()
      Go='y'
      if Go=='y':
            
            if os.path.exists(Dir+Rep)!=True: os.mkdir(Dir+Rep)
            UpFea=Dir+Rep+os.sep+Rep+'_feature_input.txt'
            UpCell=Dir+Rep+os.sep+Rep+'_feature_input_cell.txt'
            UpGene=Dir+Rep+os.sep+Rep+'_feature_input_gene.txt'
            shutil.copy2(GeneIn[:-9]+'_cell.txt',UpCell)
            print (Dir+Rep)
           
            print (Fea)
          #  print ('prune hit',len(HitPos),HitPos[:10])
            Fea1=Fea.drop(map(int,HitPos),axis=1)
            print (Fea1)
            Fea1.to_csv(UpFea,header=False,index=False)   
            GeneLs=open(GeneIn,'r').readlines()
            GeneLs1=np.array(GeneLs)
            GeneLs2=np.delete(GeneLs1,HitPos)
            OutF=open(UpGene,'w')
            OutF.write(''.join(list(GeneLs2)))
            OutF.close()
            return UpFea,UpCell,UpGene               
