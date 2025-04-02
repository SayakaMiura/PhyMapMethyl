import sys
import glob
import os
import pandas as pd
import subprocess
import numpy
import Functions
from input.Params import Params



InputType=sys.argv[2] 
Cut=float(sys.argv[11])       
PreAllHit=sys.argv[10]


python='python3'
def GetOut(OutFile,In):
    OutF=open(OutFile,'w')
    OutF.write(In)
    OutF.close()

Arg=Params()
Dir,ID,ESLin,ProInfo,Op,FeaIn0,Clo2Cell=Arg.parse(sys.argv)
if sys.argv[2]=='-Tree': 
    Op='Sample'
    ESLin=sys.argv[3][:-4]+'.txt'
    CellOr=[]
    Cell2Tu={}
elif sys.argv[2]=='-TreeC': Op='Cell'
else: open('a','r').readlines()
print (Dir,ID,ESLin,ProInfo,Op,FeaIn0) 

print ('get group info')
BP2Genels,GeneGroup,GroupLasso=Arg.GetGroupInfo(sys.argv)       
print ('Group:',GeneGroup,len(BP2Genels),GroupLasso)

LambdaArg,Rank,zLs,yLs,ID1=Arg.GetLambdaSet(sys.argv,ESLin)

print ('Response:',ESLin,zLs,yLs,Rank)
print ('ID,ID1',ID,ID1)

Pro2In,Entr2Pro,ProHea=Arg.GetGeneinInfo(ProInfo)
print ('making ESL input files')
print (Op)

ESLRes,TuOr,FeaIn,GroIn,MapIn, FeildIn=Arg.MakeESLInputFiles(ESLin,Dir,BP2Genels,GroupLasso,Entr2Pro,ID1,Op,FeaIn0,ID)
if Op=='Cell':
    CellOr,Cell2Tu=Arg.GetTuCellLs(Clo2Cell,FeaIn0[:-4]+'_cell.txt',Op)
    print ('cell count',len(CellOr),len(Cell2Tu))
OutGroupID=''
HitFeaNodeTa=Dir+os.sep+ID1+'_HitGeneNode.txt'
if os.path.exists(FeaIn)==True and os.path.exists(HitFeaNodeTa)!=True:            
   RowC=len(ESLRes)
   Ind=0
 
   AllPhyScore=['\t'.join(['ID','Response (y)','Sample','yhat','Score (y - yhat)'])+'\n']   
   Gene2ShiftNode={}
   Gene2Pos={} 

   while Ind<RowC:
    #print (Ind,RowC)
    Line=ESLRes.loc[[Ind]]

    IDn=str(list(Line.loc[:,'ID'])[0])
    Go='y'

    if IDn.strip()!='nan' and IDn.strip()!='Outgroup':
           Ingroup=list(Line.loc[:,'Set1'])[0].strip().split(';')
           Outgroup=list(Line.loc[:,'Set2'])[0].strip().split(';')
           All=Ingroup+Outgroup
           All=list(set(All))
  
                   
    else: 
       Go='n'
       if  IDn.strip()=='Outgroup':
            OutGroupID= list(Line.loc[:,'Set1'])[0].strip().split(';')[0]
            OutGroupID=OutGroupID.replace(OutGroupID.split('-')[0]+'-','')
            print ('outgroup ID',OutGroupID)  
          
    if Go=='y':
     ResIn,ResIn1,SweiIn=Arg.MakeResponseFile(CellOr,Cell2Tu,TuOr,Ingroup,Outgroup,Op,Dir,IDn,ID1)
 
     if Op=='Sample':
         ResOr,SampHead=Arg.GetInfo(ResIn,TuOr)                 
     elif Op=='Cell':
         ResOr,SampHead=Arg.GetInfo(ResIn,CellOr)  

     Rank='n'
     if os.path.exists(SweiIn)==True: #Rank=='n':
       print ('run ESL',zLs,yLs)
       for z in zLs:

         for y in yLs: 

            Out=Dir+os.sep+IDn+'_'+ID1+'_'+'out_feature_weight'+str(z)+'_'+str(y) 
            os.system(' '.join(['python3 run_SLEP.py',Out,FeaIn,SweiIn,GroIn,ResIn,MapIn]))
            if os.path.exists(Out+'.txt')==True:
                 WeiTa=Out+'.txt'
                 Gene2ShiftNode,Gene2Pos=Functions.SumSeletectGenes(WeiTa,Gene2ShiftNode,Gene2Pos,IDn,Cut)            
            else: print ('SLEP was failed')
            if os.path.exists(Out+'_temp.txt')==True: os.remove(Out+'_temp.txt')

            if os.path.exists(Out+'.xml')==True: os.remove(Out+'.xml')
              
     if os.path.exists(ResIn)==True: os.remove(ResIn)  
     if os.path.exists(ResIn[:-4]+'1.txt')==True: os.remove(ResIn[:-4]+'1.txt') 
     if os.path.exists(SweiIn)==True: os.remove(SweiIn)                           


    Ind+=1  
  

   print ('output selected genes...',len(Gene2ShiftNode))
   
   if Op=='Sample':
       BetaTa=pd.read_csv(FeaIn,sep=',',header=None)
       Functions.mergeFeatureValue(BetaTa,TuOr,Gene2ShiftNode,Gene2Pos,HitFeaNodeTa,Op)
       BetaTa=''
   elif Op=='Cell':
       Functions.mergeFeatureValueCell_simple(CellOr,Gene2ShiftNode,Gene2Pos,HitFeaNodeTa,Clo2Cell) 
   

