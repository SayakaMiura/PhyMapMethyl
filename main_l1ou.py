import sys
import glob
import os
import pandas as pd
import subprocess
import numpy
import Functions
from input.Params import Params



python='python3'
def GetOut(OutFile,In):
    OutF=open(OutFile,'w')
    OutF.write(In)
    OutF.close()
def Merge(SuppTa,PreAllHit,BooCut):
   Out=SuppTa[:-4]+'1.txt'
   OutB=SuppTa[:-4]+'2.txt'
   Ta=pd.read_csv(SuppTa,sep='\t')
   IDs=Ta['Gene']+'\t'+Ta['Position from 0'].astype(str)
   Add=[]
   PreTa=open(PreAllHit,'r').readlines()[1:]
   for i0 in PreTa:
       i=i0.split('\t')
       ID=i[0]+'\t'+i[1]
       if ID not in IDs: Add.append(i0)
   print ('number of hit',len(Add))
   out=open(SuppTa,'r').readlines()+Add
   GetOut(Out,''.join(out))   
   Ta1=pd.read_csv(Out,sep='\t') 
   Ta2 = Ta1.drop('Boo', axis=1)
   BooSco=Ta2['BooSco']
   InTot=len(BooSco)
   out2=[out[0]]

   out=out[1:]
   c=0
   while c<InTot:

       if BooSco[c]>=BooCut:
           out2+=[out[c]]
           
       c+=1
   GetOut(OutB,''.join(out2))      

ESLin=sys.argv[1] 
FeaIn0=sys.argv[2]
Op=sys.argv[3]  
PreAllHit=sys.argv[4]

if Op=='Sample':   
   Treeou=ESLin[:-4]+'.nwk'
   NodeIDin=ESLin[:-4]+'.txt'
   Rfile=sys.argv[5]
elif Op=='Cell':   
   Treeou=ESLin[:-12]+'.nwk'
   NodeIDin=ESLin[:-12]+'_ingroup.txt'
   BooScoCut=float(sys.argv[-1])
else: 
   Treeou=''
   NodeIDin='' 
if Op=='Sample':     
   TreeouTxt=Treeou[:-4]+'.txt'
   TreeTxt=pd.read_csv(TreeouTxt,sep='\t')
   OutGroupID= TreeTxt.loc[len(TreeTxt)-1,'Set1']
   OutGroupID=OutGroupID.split('-')[-1]  
else: OutGroupID=''    

if Op=='Cell':
    Functions.Do1lou(PreAllHit,Treeou,NodeIDin,Op,OutGroupID) #need
else:

    Functions.Do1lou_bulk(PreAllHit,Treeou,NodeIDin,Op,OutGroupID,Rfile) #need    

SuppTa=PreAllHit[:-4]+'_Supported.txt' 

if Op=='Cell':
    Merge(SuppTa,PreAllHit,BooScoCut)
if os.path.exists(SuppTa[:-4]+'2.txt')==True and Op=='Cell':

    print ('make box plot for each gene') 

    Com= 'Rscript --vanilla boxplot.r '+FeaIn0[:-4]+' '+ESLin[:-12]+' '+SuppTa[:-4]+'2.txt' #_ingroup.txt
     
    print (Com)              
    os.system(Com)

