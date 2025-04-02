from input.Loop import Loop
import shutil
import pandas as pd
import os
import sys
import glob
from input.Params import Params
import Functions
import numpy
import math

Input1=sys.argv[1]
InputType=sys.argv[2] 
InputT=sys.argv[3]
Input=InputT
InputOri=InputT
WeiCut=float(sys.argv[4])
RepMax=float(sys.argv[5]) 
MaxSigleHit=float(sys.argv[6]) 
MinSigleHit=float(sys.argv[7]) 
GeneInfo=sys.argv[8]
Range=['-Range','0.2','0.2','0','0']
Min=float(Range[1])
Max=float(Range[2])
LamLs=list(numpy.logspace(math.log10(Min),math.log10(1), 10, endpoint=True))

Out=Input[:-4]+'_DetectedGenes.txt'

if Input1[-4:]!='.txt':
    OriDir=Input1+os.sep
    Op='Sample'
    Rfile='Runl1ouWObooBulk.r'
    if GeneInfo=='RNA':  Rfile='Runl1ouWObooBulkRNA.r'
    print ('making feature input files')
    Arg=Params()
    ID1='input'
    ESLRes,TuOr,FeaIn,GroIn,MapIn, FeildIn=Arg.MakeESLInputFiles_simple(InputT[:-4]+'.txt',Input1,ID1)
    TotSiteIn=int(open(MapIn,'r').readlines()[-1].split('\t')[0])
    print (TuOr)

    
else:
    Op='Cell'
    OriDir=Input1.replace(Input1.split(os.sep)[-1],'')
    FeaIn=Input1
    BooScoCut=50
    MapIn=FeaIn[:-4]+'_gene.txt'
    TotSiteIn=len(open(MapIn,'r').readlines())
   
print ('total features in input',TotSiteIn)  

RunESLsrg=[FeaIn,InputType,InputT,GeneInfo]+Range
RunESLsrgini=[FeaIn,InputType,InputT,GeneInfo]+Range

Tot=0
FeaInOri=FeaIn
OriGeneIn=FeaIn[:-4]+'_gene.txt'
BooCut=60
ProcessRep=Loop()
Rep=1
AllHitNode=FeaIn.replace(FeaIn.split(os.sep)[-1],'')+Input.split(os.sep)[-1][:-4]+'_AllHitNode.txt'
PreHit=[]
PreHitFile=AllHitNode.replace(AllHitNode.split(os.sep)[-1],'')+'PhyMap_All.txt'
print (PreHitFile)
if os.path.exists(PreHitFile)==True:
        Pdf=pd.read_csv(PreHitFile,sep='\t')
        PreHit=Pdf['Gene'].to_list()
        PreHit=list(set(PreHit))
LamP=1
Repeat='y'
Comment=''
Note=''
while Min<1.0: 
 # print (MapIn)

  Go='y'
  os.system('python3 main_SLEP.py '+' '.join(RunESLsrg)+' '+AllHitNode+' '+str(WeiCut)) #AllHit is not used
  print ('python3 main_SLEP.py '+' '.join(RunESLsrg)+' '+AllHitNode+' '+str(WeiCut))
 # BetaTa=pd.read_csv(FeaIn,sep=',',header=None)
  if InputType=='-TreeC':
      ESLin=Input[:-4]+'_ingroup.txt'
  else: ESLin=Input    
  HitFeaNodeTa=FeaIn.replace(FeaIn.split(os.sep)[-1],'')+ESLin.split(os.sep)[-1][:-4]+'_HitGeneNode.txt'
  if os.path.exists(HitFeaNodeTa)!=True:
    print ('SLEP was failed. Use next lambda')
    if LamP>=len(LamLs): 
              LamLs=list(numpy.logspace(math.log10(LamLs[-1]),math.log10(LamLs[-1]+0.4), 10, endpoint=True)
)
              LamP=0     
                
    Min=LamLs[LamP]
    Max=LamLs[LamP]

    print (Min)
    RepMax+=1
    Comment+='SLEP was failed and SLEP lambda was updated: '+str(Min)+'\n'

    RunESLsrg[4]='-Range '+str(Min)+' '+str(Max)+' 0 0'
    print (RunESLsrg)
    LamP+=1       

    
  else:
    HitTa0=pd.read_csv(HitFeaNodeTa,sep='\t')
    SingleHitC=len(HitTa0['Gene'])
    print ('Rep ',Rep,', Number of features seletected by SLEP: ',SingleHitC)

    if Rep==1:     
      if SingleHitC<MaxSigleHit:    
        AllHitNodeIn=pd.read_csv(HitFeaNodeTa,sep='\t')
        SingleHitC=len(AllHitNodeIn['Gene'])
        Comment+='Initial SLEP selection: '+str(SingleHitC)+'\n'
      else:
        os.remove(HitFeaNodeTa)
        Go='n'
        if LamP>=len(LamLs): 
              LamLs=list(numpy.logspace(math.log10(LamLs[-1]),math.log10(LamLs[-1]+0.4), 10, endpoint=True)
)
              LamP=0     
                
        Min=LamLs[LamP]
        Max=LamLs[LamP]

        print (Min)
    
        Comment+='SLEP was failed and SLEP lambda was updated: '+str(Min)+'\n'

        RunESLsrg[4]='-Range '+str(Min)+' '+str(Max)+' 0 0'
        print (RunESLsrg)
        LamP+=1         

        Comment+='Initial SLEP selection was too many and lambda was updated: '+str(SingleHitC)+' '+str(Min)+'\n'        

    else:   
        AddN=pd.read_csv(HitFeaNodeTa,sep='\t')
        SingleHitC=len(AddN['Gene'])#.unique())
        if SingleHitC<MaxSigleHit:
            AllHitNodeIn = pd.concat([AllHitNodeIn, AddN])
            AllHitNodeIn = AllHitNodeIn.drop_duplicates()
            Comment+=str(Rep)+' SLEP selection: '+str(SingleHitC)+'\n'
    if Go=='y':        
     AllHitNodeIn.to_csv(AllHitNode,sep='\t',index=False)    
     GeneIn=FeaIn[:-4]+'_gene.txt'
     Dir=FeaIn.replace(FeaIn.split(os.sep)[-1],'')
     Tot+=SingleHitC
     if TotSiteIn==Tot:
       Comment+=str(Rep)+' SLEP selected all features, so the computation was terminated: '+str(Tot)+' '+str(TotSiteIn)+'\n' 
       Repeat='n' 
     if SingleHitC>MaxSigleHit: 
        RepMax+=1
        Comment+=str(Rep)+' SLEP selection is too many (discarded): '+str(SingleHitC)+'\n' 
        print ('too many hits',SingleHitC)  
        BetaTa=pd.read_csv(FeaIn,sep=',',header=None)  
        FeaIn,UpCell,UpGene =ProcessRep.UpInput_copy(BetaTa,GeneIn,HitFeaNodeTa,BooCut,'Rep'+str(Rep))  
        BetaTa='s'
        if LamP>=len(LamLs) and SingleHitC>MaxSigleHit: 
              LamLs=list(numpy.logspace(math.log10(LamLs[-1]),math.log10(LamLs[-1]+0.4), 10, endpoint=True)
)
              LamP=1         
        Min=LamLs[LamP]
        Max=LamLs[LamP]
        Comment+='SLEP lambda was updated: '+str(Min)+'\n'
        LamP+=1 
     else:
        BetaTa=pd.read_csv(FeaIn,sep=',',header=None)
        FeaIn,UpCell,UpGene =ProcessRep.UpInput0(BetaTa,GeneIn,HitFeaNodeTa,BooCut,'Rep'+str(Rep)) 
        BetaTa=''
     if SingleHitC<MinSigleHit: 
        if LamP>=len(LamLs): 
              LamLs=list(numpy.logspace(math.log10(LamLs[-1]),math.log10(LamLs[-1]+0.4), 10, endpoint=True)
)
              LamP=1
              Min=LamLs[LamP]
              Max=LamLs[LamP]  
              print ('repeat again',LamP,LamLs,Min)  
              LamP+=1                     
        else:
            Min=LamLs[LamP]
            Max=LamLs[LamP]
            LamP+=1  
        Comment+='SLEP lambda was updated: '+str(Min)+'\n' 
     if SingleHitC==0: RepMax+=1       
     if  Rep>RepMax:

          CHit=AllHitNodeIn['Gene'].to_list()
          AllHitR=PreHit+CHit
          AllHitU=list(set(AllHitR))
          if len(AllHitU)==len(PreHit): RepMax+=1         
          else:Repeat='n'
     if Repeat=='y':#FeaIn!='Done':
        Dir=UpGene.replace(UpGene.split(os.sep)[-1],'')
        UpInput=Dir+Input.split(os.sep)[-1]
        shutil.copy2(Input,UpInput)
        shutil.copy2(Input[:-4]+'.txt',UpInput[:-4]+'.txt')       
        if Op=='Sample': shutil.copy2(OriDir+'input_feature_input_sample.txt',Dir+'Rep'+str(Rep)+'_feature_input_sample.txt')                           
        RunESLsrg=[FeaIn,InputType,UpInput,GeneInfo,'-Range '+str(Min)+' '+str(Max)+' 0 0']
        Input=UpInput
        Rep+=1
 
if Tot==0:
   Comment+='There is no shift detected by SLEP.\n' 
else:
    print ('make l1ou input')

    OriGeneOrder=open(OriGeneIn,'r').readlines()
    AllHitNodeIn=open(AllHitNode,'r').readlines()
    out=[AllHitNodeIn[0]]
    AllHitNodeIn=AllHitNodeIn[1:]
    Gene2NodeLs={}
    Gene2beta={}
    for i in AllHitNodeIn:
       i=i.split('\t')
       Gene=i[0]
       if Gene not in PreHit:
        NodeLs=i[2].split(';')
        Gene2NodeLs[Gene]=Gene2NodeLs.get(Gene,[])+NodeLs
        Gene2beta[Gene]=i[3:]
    for Gene in Gene2NodeLs:
        NodeLs=list(set(Gene2NodeLs[Gene]))
        NodeLs=';'.join(NodeLs)
        Pos=OriGeneOrder.index(Gene+'\n')
        AllIn=[Gene,str(Pos),NodeLs]+Gene2beta[Gene]
        out+=['\t'.join(AllIn)]    
    OutF=open(AllHitNode[:-4]+'1.txt','w')
    OutF.write(''.join(out))
    OutF.close()    
    if InputType=='-TreeC': pass     
    else: 
        print ('run l1ou')
        Com='python3 main_l1ou.py '+InputOri+' NA Sample '+AllHitNode[:-4]+'1.txt '+Rfile
        os.system(Com) 
        print (Com) 
        if os.path.exists(AllHitNode[:-4]+'1_Supported.txt')==True:
            OutF,Note=Functions.ComputeFP(AllHitNode[:-4]+'1.txt',AllHitNode[:-4]+'1_Supported.txt')
        else: 
            Comment+='There is no shift detected by l1ou.\n' 
            Note='' 
            OutF=AllHitNode[:-4]+'1_Note.txt'     

if Go=='y':
  NoteIn=['Final SLEP lambda: ',str(WeiCut),'SLEP maximum selection: ',str(MaxSigleHit),'SLEP minimum selection: ',str(MinSigleHit),'SLEP iteration: ',str(RepMax)]  
  Note+='\n'.join(NoteIn)+'\n\n'+Comment
  Functions.GetOut(OutF,Note)                
   
     
        
        
        
    
            
