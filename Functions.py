import sys
import glob
import os
import pandas as pd
import subprocess
import numpy
from random import choices

def GetGeneinInfo(ProInfo):
 Pro2In={}
 Entr2Pro={}
 ProHea=''
 if os.path.exists(ProInfo)==True:
    print ('reading probe information..')
    df=pd.read_csv(ProInfo,sep='\t')
    print (df)
    ProInfo=open(ProInfo,'r').readlines()
    
    ProHea='\t'+ProInfo[0].strip()
    ProInfo=ProInfo[1:]
    IDls=[]
    for i in ProInfo:
        
        Pro2In[i.split('\t')[0]]='\t'+i.strip()
        IDls.append(i.split('\t')[0])
    Entre=list(df.loc[:,'entrezgene'])
    c=0
    while c<len(IDls):
        if str(Entre[c])!='nan' and str(Entre[c]).find(';')==-1:
            Entr2Pro[str(int(Entre[c]))]=IDls[c]
        c+=1
    print (len(Entr2Pro))
    #open('a','r').readlines() :  
 return Pro2In,Entr2Pro,ProHea 
def GetTuOrder(ESLRes):
    Inres0=list(ESLRes.loc[:,'Set1'])
    Outres0=list(ESLRes.loc[:,'Set2'])
    Allres0=Inres0+Outres0
    TuOr=[]
    for i in Allres0:
      if str(i)!='nan':
        i=i.split(';')
        TuOr+=i
    TuOr=list(set(TuOr))
    print (len(TuOr))
    return TuOr 
def MakeFeatureFilePreprocess(Dir,TuOr):
    s=0
    Amat=[]
    FileLs=[]
    for Tu in TuOr:
       File=Dir+os.sep+Tu+'.txt'
       if os.path.exists(File)!=True: 
           print ('sample feature file does not exist',File)
       else:
           
           Beta=pd.read_csv(File,sep='\t',skiprows=[0,1])
           if s!=0:
               if ProbLs!=list(Beta.loc[:,'ID_REF']):
                   print ('the order of ID_REF is different, please fix')
                   break
               else: FileLs.append(File)
           else: FileLs.append(File)        
           s+=1        
           ProbLs=list(Beta.loc[:,'ID_REF'])
           Bls0=list(Beta.loc[:,'VALUE'])
           Bls=[]
           for B in Bls0:
              # if B<Cut: Bls.append(0)
              # else: 
              Bls.append(B)
           Amat.append(Bls)#(','.join(map(str,Bls))+'\n')

    return Amat,FileLs, ProbLs  
def MakeMapIn(ProbNum,ProbLs1,MapIn):     
       Map=['0\t\n']
       c=0
       while c<ProbNum:
           Map.append(str(c+1)+'\t'+ProbLs1[c].strip()+'\n')
           c+=1
       GetOut(MapIn,''.join(Map).replace(' ','')) 
def GetBP2genels(BP2Genels,Entr2Pro,ProbLs):
           BP2genels={}
           for BP in BP2Genels:
               Genels=BP2Genels[BP]

           
               for G in Genels:
                   Pro=Entr2Pro.get(G,'')
               
                   if Pro in ProbLs:

                       BP2genels[BP]=BP2genels.get(BP,[])+[Pro]   
           return BP2genels                                
#def MakeESLInputFiles(ESLin,Dir,BP2Genels,GroupLasso,Entr2Pro,ID1):
#  ESLRes=pd.read_csv(ESLin,sep='\t')
#  print (ESLRes)
#  TuOr=GetTuOrder(ESLRes)
##    Inres0=list(ESLRes.loc[:,'Set1'])
##    Outres0=list(ESLRes.loc[:,'Set2'])
##    Allres0=Inres0+Outres0
##    TuOr=[]
##    for i in Allres0:
##      if str(i)!='nan':
##        i=i.split(';')
##        TuOr+=i
##    TuOr=list(set(TuOr))
##    print (len(TuOr))
#    ####
##open('a','r').readlines()
#  FeaIn=Dir+os.sep+ID1+'_feature_input.txt'
#  MapIn=Dir+os.sep+ID1+'_feature_mapping_input.txt' 
#  GroIn=Dir+os.sep+ID1+'_group_indices_input.txt'    
#  FeildIn=Dir+os.sep+ID1+'_field_input.txt'
#  if os.path.exists(FeaIn)!=True:
#    Amat,FileLs,ProbLs=MakeFeatureFilePreprocess(Dir,TuOr)
##    s=0
##    Amat=[]
##    FileLs=[]
##    for Tu in TuOr:
##       File=Dir+os.sep+Tu+'.txt'
##       if os.path.exists(File)!=True: 
##           print ('sample feature file does not exist',File)
##       else:
#           
##           Beta=pd.read_csv(File,sep='\t',skiprows=[0,1])
##           if s!=0:
##               if ProbLs!=list(Beta.loc[:,'ID_REF']):
##                   print ('the order of ID_REF is different, please fix')
##                   break
##               else: FileLs.append(File)
##           else: FileLs.append(File)        
##           s+=1        
##           ProbLs=list(Beta.loc[:,'ID_REF'])
##           Bls0=list(Beta.loc[:,'VALUE'])
##           Bls=[]
##           for B in Bls0:
##              # if B<Cut: Bls.append(0)
##              # else: 
##              Bls.append(B)
##           Amat.append(Bls)#(','.join(map(str,Bls))+'\n')
#    PathIn=Dir+os.sep+'path.txt'
#    GetOut(PathIn,'.fas\n'.join(ProbLs)+'.fas\n') 
#    
#    if len(TuOr)==len(FileLs):
#      # FeaIn=Dir+os.sep+ID1+'_feature_input.txt'
#       print ('remove probes with nan. Making table of beta values...')
#       FeaTa=pd.DataFrame(Amat,columns=ProbLs)#pd.read_csv(FeaIn0,sep=',',names=ProbLs)
#       print ('pruning')
#       FeaTaClean=FeaTa.dropna(axis=1)
#       FeaTaClean.to_csv(FeaIn,sep=',',header=False,index=False)
#       ProbLs1=list(FeaTaClean.columns.values)
#       ProbNum=len(ProbLs1)   
#       print (ProbNum,len(ProbLs),ProbLs1[0])
#    #   MapIn=Dir+os.sep+ID1+'_feature_mapping_input.txt' 
#       MakeMapIn(ProbNum,ProbLs1,MapIn)     
#  #     Map=['0\t\n']
#  #     c=0
#  #     while c<ProbNum:
#  #         Map.append(str(c+1)+'\t'+ProbLs1[c].strip()+'\n')
#  #         c+=1
#  #     GetOut(MapIn,''.join(Map).replace(' ',''))
#       #####
#    #   GroIn=Dir+os.sep+ID1+'_group_indices_input.txt'    
#    #   FeildIn=Dir+os.sep+ID1+'_field_input.txt'

       
#       PosLs=list(range(1,ProbNum+1))
#       NumStr=','.join(map(str,PosLs))
#       Seq=list(range(1,ProbNum+1))
#       Seq1=list(map(str,Seq))
#       PosAllSep=','.join(Seq1)
#       AllOne=['1']*ProbNum
#       AllOne1=','.join(AllOne)
#       print (len(AllOne),ProbNum,len(Seq),len(Seq1),PosAllSep.count(','),AllOne1.count(','))
#       if BP2Genels!={}:
#           BP2genels=GetBP2genels(BP2Genels,Entr2Pro,ProbLs)
#       #    BP2genels={}
#       #    for BP in BP2Genels:
#       #        Genels=BP2Genels[BP]

           
#       #        for G in Genels:
#       #            Pro=Entr2Pro.get(G,'')
#               
#       #            if Pro in ProbLs:

#       #                BP2genels[BP]=BP2genels.get(BP,[])+[Pro]   
###############                       
#       if GroupLasso=='n':#BP2Genels=={}: 
#           print (GroIn)
#         #  GetOut(GroIn,PosAllSep+'\n'+PosAllSep+'\n'+AllOne1+'\n')
#           GetOut(GroIn,'1\n'+str(ProbNum)+'\n'+str(ProbNum**0.5)+'\n')  
#      # open('a','r').readlines() 
      
#       else:  
#           Gdf={}
#           ProPosOrAll=[]
#           Minor=[]
#           BP2genels={}
#           AllBPgeneLs=[]
#           S=1
#           for BP in BP2Genels:
#               Genels=BP2Genels[BP]
#               Genels1=[]
#               ProPosOr=[]
#           
#               for G in Genels:
#                   Pro=Entr2Pro.get(G,'')
#               
#                   if Pro in ProbLs:
#               #    print (Pro,ProbLs.index(Pro))
#                       ProPosOr.append(ProbLs.index(Pro)+1)
#                       BP2genels[BP]=BP2genels.get(BP,[])+[Pro]
#                       AllBPgeneLs.append(Pro)
#               if len(ProPosOr)<2: 
#                  Minor+=ProPosOr
#                  print ('<2',BP,Genels)
#              #open('a','r').readlines()
#               else: 
#                   ProPosOrAll+=ProPosOr
#                   E=S+len(ProPosOr)-1
#                   Gdf[BP]=[str(S),str(E),(float(len(ProPosOr)))**0.5]
#                   S=E+1    
#                #   open('a','r').readlines()
#           Miss=[]
#           c=1
#           while c<=len(ProbLs):
#               if c not in ProPosOrAll: Miss.append(c)
#               c+=1
#           print (len(Miss),c)
#       #open('a','r').readlines()    
#           ProPosOrAll+=Miss    
#           E=S+len(Miss)-1    
#           Gdf['Miss']=[str(S),str(E),(float(len(Miss)))**0.5] 
#           GetOut(FeildIn,','.join(map(str,map(int,ProPosOrAll))))
#           Gdf1=pd.DataFrame(Gdf)
#           print (Gdf1)
#           Gdf1.to_csv(GroIn[:-4]+'_withHeader.csv',index=False)
#           Gdf1.to_csv(GroIn,index=False,header=False)    
#  return ESLRes,TuOr,FeaIn,GroIn,MapIn, FeildIn          
def SumSeletectGenes(WeiTa,Dic,DicP,Node,Cut):
    if len(open(WeiTa,'r').readlines())>0:
     PosMethWei=pd.read_csv(WeiTa, header=None,sep='\t')
     print (PosMethWei)
     PosMethWei.columns =['PositionFrom1','Probe','Weight'] 
     HitLs=PosMethWei.loc[:,'Probe']
     TarPos=PosMethWei.loc[:,'PositionFrom1']  
     Wei=PosMethWei.loc[:,'Weight'] 
     c=0  
     for Hit in HitLs:
      if isinstance(Wei[c], int)==True or isinstance(Wei[c], float)==True:
           if abs(Wei[c])>Cut:
              Dic[Hit]=Dic.get(Hit,[])+[Node]
              DicP[Hit]=int(TarPos[c])-1
      c+=1
    print ('Weight cutof, the number of detected shifts',Cut,len(DicP))  
    return Dic,DicP    
def ComputePhyScore(WeiTa,BetaTa,SampHead): #WeiTa[:-4]+'_PhyScore.txt'
            
            

                 PosMethWei=pd.read_csv(WeiTa, header=None,sep='\t')
                 PosMethWei.columns =['PositionFrom1','Probe','Weight'] 
                 print (PosMethWei)
                 TarPos=PosMethWei.loc[:,'PositionFrom1']
                 print (TarPos[0],len(TarPos))
 
                 print ('extracting hit probes beta values...')
                 c=0
                 Samp2SumY={}
                 for Pos in TarPos:
                     Pos=int(Pos)-1
                     MethID=PosMethWei.loc[:,'Probe'][c]
                     WeiV=PosMethWei.loc[:,'Weight'][c]
                     BetaLs=list(BetaTa.iloc[:,Pos])
                     T=0
                     while T<len(SampHead):
                         Samp2SumY[SampHead[T]]=Samp2SumY.get(SampHead[T],0)+(BetaLs[T]*WeiV)
                         T+=1
                     c+=1
                 print ('making output...')
                 print ('compute phyScore...')
                 out=['\t'.join(['ID','Response (y)','Sample','yhat','Score (y - yhat)'])+'\n']
                 for Samp in Samp2SumY:
                     SumY=Samp2SumY[Samp]
                     Res=float(Samp.split('_')[0])
                     Score=Res-SumY
                     out.append('\t'.join([WeiTa.split(os.sep)[-1],str(Res),Samp[(len(str(Res))+1):],str(SumY),str(Score)])+'\n')
 
                 return out #['Response (y)','Sample','yhat','Score (y - yhat)']
def  GetInfo(ResIn,TuOr):                 
     ResOr=open(ResIn,'r').readlines()
     SampHead=[]
     c=0
     while c<len(ResOr):
                     SampHead.append(ResOr[c].strip()+'_'+TuOr[c])
                     c+=1
     print (SampHead) 
       
     return ResOr,SampHead      
def mergeFeatureValue(BetaTa,SampHead,GeneDic,Gene2Pos,OfileName,Op): #Gene2Pos from 0
          SampHead0=[]
          for i in SampHead:
               if Op=='Sample': SampHead0.append(i[(len(i.split('-')[0])+1):])
               elif Op=='Cell': SampHead0.append(i)
          out=['Gene\tPosition from 0\tNode\t'+'\t'.join(SampHead0)+'\n']#'\t'.join(ProbLsIn)+'\t'+'\n']  
          for Gene in GeneDic:
              NodeLs=list(set(GeneDic[Gene]))
              NodeLs.sort()
              Pos0=Gene2Pos[Gene]
              BetaLs=list(BetaTa.iloc[:,Pos0])
              out.append(str(Gene)+'\t'+str(Pos0)+'\t'+';'.join(NodeLs)+'\t'+'\t'.join(list(map(str,BetaLs)))+'\n')
 
          
          GetOut(OfileName,''.join(out))  
def mergeFeatureValueCell(BetaTa,SampHead,GeneDic,Gene2Pos,OfileName,Clo2Cell): #Gene2Pos from 0
          SampHead0=list(Clo2Cell.keys())
      #    for i in SampHead:
      #         if Op=='Sample': SampHead0.append(i[(len(i.split('-')[0])+1):])
      #         elif Op=='Cell': SampHead0.append(i)
          out=['Gene\tPosition from 0\tNode\t'+'\t'.join(SampHead0)+'\n']#'\t'.join(ProbLsIn)+'\t'+'\n'] 
          outBoo=['Gene\tPosition from 0\tNode\t'+'\t'.join(SampHead0)+'\n']
          print (len(GeneDic),SampHead0) 
          for Gene in GeneDic:
            #  print (Gene)
              NodeLs=list(set(GeneDic[Gene]))
              NodeLs.sort()
              Pos0=Gene2Pos[Gene]
              BetaLs=list(BetaTa.iloc[:,Pos0])
              outIn=[str(Gene),str(Pos0),';'.join(NodeLs)]
          
              for Clone in SampHead0:
                  CellLs=Clo2Cell[Clone]
                #  if Gene=='Clu': print (len(CellLs))
                  BetaClo=[]
                  for Cell in CellLs:
                      Ind=SampHead.index(Cell)
                      BetaClo.append(float(BetaLs[Ind]))
                     # if Gene=='Clu': print (Cell,float(BetaLs[Ind]))
                #  if Gene=='Clu': print (BetaClo)    
                #  print (Clone,CellLs,BetaClo)
                  Median=numpy.median(BetaClo) #numpy.mean(BetaClo)#
                  
                 # print (Median)
                  outIn.append(str(Median))    
              out.append('\t'.join(outIn)+'\n')
              BooNum=100
              BooID=0
              while BooID<BooNum:
               outInBoo=[str(Gene),str(Pos0),';'.join(NodeLs)]
          
               for Clone in SampHead0:
                  CellLs=Clo2Cell[Clone]
                #  if Gene=='Clu': print (len(CellLs))
                  BetaClo=[]
                  for Cell in CellLs:
                      Ind=SampHead.index(Cell)
                      BetaClo.append(float(BetaLs[Ind]))
                     # if Gene=='Clu': print (Cell,float(BetaLs[Ind]))
                #  if Gene=='Clu': print (BetaClo)    
                #  print (Clone,CellLs,BetaClo)
                  BooSamp=choices(BetaClo, k=len(BetaClo))
                  Median=numpy.median(BooSamp) #numpy.mean(BetaClo)#
                  
                 # print (Median)
                  outInBoo.append(str(Median))    
               outBoo.append('\t'.join(outInBoo)+'\n')              
               BooID+=1    
         # open('a','r').readlines()       
          GetOut(OfileName,''.join(out))   
          GetOut(OfileName[:-4]+'_boo.txt',''.join(outBoo))  
def mergeFeatureValueCell_simple(SampHead,GeneDic,Gene2Pos,OfileName,Clo2Cell): #Gene2Pos from 0
          SampHead0=list(Clo2Cell.keys())

          out=['Gene\tPosition from 0\tNode\t'+'\t'.join(SampHead0)+'\n']#'\t'.join(ProbLsIn)+'\t'+'\n'] 

          print (len(GeneDic),SampHead0) 
          for Gene in GeneDic:

              NodeLs=list(set(GeneDic[Gene]))
              NodeLs.sort()
              Pos0=Gene2Pos[Gene]

              outIn=[str(Gene),str(Pos0),';'.join(NodeLs)]
          

              out.append('\t'.join(outIn)+'\n')
     
          GetOut(OfileName,''.join(out))   
     
def makel1ouin(InFile,BetaTa,SampHead,Clo2Cell):
          OfileName=InFile[:-4]+'1.txt'
          AllHit=open(InFile,'r').readlines()
          SampHead0=AllHit[0].strip().split('\t')[3:]#list(Clo2Cell.keys())
          print (SampHead0)
          outBoo=[AllHit[0]]
          out=[AllHit[0]]
          AllHit=AllHit[1:]
         # out=['Gene\tPosition from 0\tNode\t'+'\t'.join(SampHead0)+'\n']#'\t'.join(ProbLsIn)+'\t'+'\n'] 
        #  outBoo=['Gene\tPosition from 0\tNode\t'+'\t'.join(SampHead0)+'\n']
          print (len(AllHit),SampHead0) 
          print ('generating bootstrap replicates')
          for i in AllHit:
              Pos0=int(i.split('\t')[1])
         # for Gene in GeneDic:
            #  print (Gene)
            #  NodeLs=list(set(GeneDic[Gene]))
            #  NodeLs.sort()
            #  Pos0=Gene2Pos[Gene]
              BetaLs=list(BetaTa.iloc[:,Pos0])
              outIn=[i.strip()]#[str(Gene),str(Pos0),';'.join(NodeLs)]              
          
              for Clone in SampHead0:
                  CellLs=Clo2Cell[Clone]
                  BetaClo=[]
                  for Cell in CellLs:
                    if Cell in SampHead:
                      Ind=SampHead.index(Cell)
                      BetaClo.append(float(BetaLs[Ind]))
                  Median0=numpy.median(BetaClo) #numpy.mean(BetaClo)#
                  outIn.append(str(Median0))  
              out.append('\t'.join(outIn)+'\n')  
              BooNum=100
              BooID=0

              while BooID<BooNum:
               outInBoo=[i.strip()]
          
               for Clone in SampHead0:
                  CellLs=Clo2Cell[Clone]
                  BetaClo=[]
                  for Cell in CellLs:
                     if Cell in SampHead:
                      Ind=SampHead.index(Cell)
                      BetaClo.append(float(BetaLs[Ind]))
                  BooSamp=choices(BetaClo, k=len(BetaClo))
                  Median=numpy.median(BooSamp) #numpy.mean(BetaClo)#
                  outInBoo.append(str(Median))  
               #   print (numpy.median(BetaClo),numpy.median(BooSamp))
                  
  
                                 
               outBoo.append('\t'.join(outInBoo)+'\n')              
               BooID+=1       
  
          GetOut(OfileName,''.join(out))  
          print ('done') 
          GetOut(OfileName[:-4]+'_boo.txt',''.join(outBoo))    
        #  open('a','r').readlines()                  
def generateBooInput(In,Boo):
    Input=open(In,'r').readlines()
    outBoo=[Input[0]]
    Input=Input[1:]
    for i in Input:
        outBoo.append(i*Boo)
    GetOut(In[:-4]+'_boo.txt',''.join(outBoo))                   
def AddGeneInfo(Ta,Pro2In,ProHea):
    OutFile=Ta[:-4]+'_gene.txt'
    Ta=open(Ta,'r').readlines()
    out=[Ta[0].strip()+'\t'+ProHea+'\n']
    Ta=Ta[1:]
    for i in Ta:
        Gene=i.split('\t')[0]
        out.append(i.strip()+'\t'+Pro2In.get(Gene,'')+'\n')
    GetOut(OutFile,''.join(out))          
def GetNodeID2TipLs(Ta,Gene):
   # print (Ta)
   # print (Gene)
    SubTa=Ta[Ta["rep.Hit.Gene..1...length.eModel.tree.tip.label.."] == Gene]
   # print (SubTa)
    NodeIDls=list(SubTa)
    Node2TipLs={}
    TipName=SubTa['eModel.tree.tip.label'].tolist()
    for N in NodeIDls:
        
        if N!='rep.Hit.Gene..1...length.eModel.tree.tip.label..' and N!='eModel.tree.tip.label':
         #   print (N)
            TipCode=SubTa[N].tolist()
            
          #  print (TipCode)
            Len=len(TipCode)
            TipLs=[]
            c=0
            while c<Len:
                if int(TipCode[c])==1: TipLs.append(TipName[c])
                c+=1
           # print (TipLs)
            TipLs=list(set(TipLs))
            TipLs.sort()
            Node2TipLs[N]=TipLs
   # print (Node2TipLs)
    if Node2TipLs=={}:
        Node2TipLs['NA']=['NA']
    return Node2TipLs  
def AddColumn(Ori0,AddDic):#
                    
                 #   print (Ori0)
                    Ori={}
                    for O0 in Ori0:
                        In000=Ori0[O0]
                      #  print (O0,In000)
                        #for O0ind in O0:
                        Ori[O0]=[In000[list(In000.keys())[0]]]
                    for i in AddDic:    
                        Ori[i]=AddDic[i]
                 #   print (Ori)
                    Add=pd.DataFrame(Ori)
                    return Add    
def CountHit(File):
    File=open(File,'r').readlines()[1:]
    ID2C={}
    BraLs=[]
    for i in File:
        i=i.split('\t')
        ID=i[0].strip()+'_'+i[2].strip()
        ID2C[ID]=ID2C.get(ID,0)+1
        BraLs.append(i[0].strip())
    BraLs=list(set(BraLs))    

    return ID2C,BraLs    
def MergeHit(ID2C,File):
    In=open(File,'r').readlines()
  #  out=['Gene\tPosition from 0\tNode\tBoo\n']
    out=[In[0].strip()+'\tBoo\n']
    In=In[1:]
    for i0 in In:
        i=i0.strip().split('\t')
        ID=i[0].strip()+'_'+i[2].strip()    
        C=ID2C.get(ID,0)
        out.append(i0.strip()+'\t'+str(C)+'\n')
    OutF=open(File,'w')
    OutF.write(''.join(out))
    OutF.close()   
 
def Do1lou_bulk(HitFeaNodeTa,Nwk,NodeIDin,Op,OutGroup,Rfile):

    if OutGroup!='':
        Tre=open(Nwk,'r').readlines()[0].strip()
        Tre='('+Tre.replace(';','')+':1,('+OutGroup+':0,N1:0):1);'
        Nwk=Nwk[:-4]+'1.nwk'
        GetOut(Nwk,Tre)

    if os.path.exists(HitFeaNodeTa+'_Shift.txt')!=True:
        print ('Rscript --vanilla '+Rfile+' '+HitFeaNodeTa+' '+Nwk)

        os.system('Rscript --vanilla '+Rfile+' '+HitFeaNodeTa+' '+Nwk)
        


    print ('summarize')

    Hit=HitFeaNodeTa+'_Shift.txt'
    NodeID=HitFeaNodeTa+'_TreeNode.txt'
    NodeID=pd.read_csv(NodeID,sep='\t')
    Hit=pd.read_csv(Hit,sep='\t')
    In=pd.read_csv(HitFeaNodeTa,sep='\t')
    In.insert(2, "AlphaShift", ['NA']*len(In['Node']), True)

    MaxSamC=len(In.columns)-3-2

    NodeIDin=open(NodeIDin,'r').readlines()[1:]
    NodeID2TipLsIn={}
    for i in NodeIDin:
      if i.strip()!='':
        i=i.split('\t')

        ID=i[0]
        Ls=i[1].split(';')
        In0=[]
        for ii in Ls:
            if Op=='Sample':ii=ii[(len(ii.split('-')[0])+1):]
            In0.append(ii)
        In0.sort()    
        NodeID2TipLsIn[ID]=In0 

               
    c=0

    RowN=len(In)
    outS='y'
    OutUnS='y'
    print (In)
    print (Hit)
    print ('summarizing...')
    while c<RowN:

        Gene=In["Gene"].tolist()[c]
        
        NodeIDInLs=In["Node"].tolist()[c].split(';')
        NodeID2TipLs=GetNodeID2TipLs(NodeID,Gene)
        HitBra=Hit[Hit['Gene']==Gene]

        if len(HitBra)>0:

            HitBraLs=HitBra['branch'].tolist()

            Entered='n'
            BadHit='n'
            Hitc=0
            for Bra in HitBraLs:

                TipLs=NodeID2TipLs['X'+str(Bra)]

                if len(TipLs)>1 and MaxSamC>=len(TipLs):
                        AS=HitBra["S.v"].tolist()[Hitc]
                        Ori0=In[In["Gene"]==Gene]#.to_dict()
                        Add=Ori0.head(1).copy()
                        Add.at[c, 'Node'] = ''.join(TipLs)   
                        Add.at[c, "AlphaShift"]= AS     
                        if outS=='y':
                            Entered='y'
                            outS='n'  
                            out1=Add
                        else:
                            Entered='y'
                            out1=pd.concat([out1,Add])
                Hitc+=1            
        c+=1
    if outS=='y': print ('l1ou did not detect any shifts')
    else:
        out1.to_csv(HitFeaNodeTa[:-4]+'_Supported.txt',sep='\t', index=False)   


def ComputeFP(AllHitFile,SuppFile):
    All=pd.read_csv(AllHitFile,sep='\t')  
    AllH=All['Gene'].nunique()
    Supp=pd.read_csv(SuppFile,sep='\t')  
    SuppH=Supp['Gene'].nunique()
    FP=AllH-SuppH
    FPr=1.0*FP/AllH
    print ('all,supp,FP,FP rate: ',AllH,SuppH,FP,FPr)
    out='All hit: '+str(AllH)+'\nSupported hit: '+str(SuppH)+'\nFP rate: '+str(FPr)+'\n'
    GetOut(AllHitFile[:-4]+'_Note.txt',out)
    print (out)
    return AllHitFile[:-4]+'_Note.txt',out                   
def AdjustSupport(FileS):
    Out=FileS[:-5]+'.txt'
    FileU=FileS[:-14]+'.txt'
    U=pd.read_csv(FileU,sep='\t')   
    IndexLs=list(U.index)
    out=open(FileS,'r').readlines()
    out.append('\n')
    for i in IndexLs:
        ExpreLs=U.iloc[i][3:]
        HitLs=U.iloc[i][2].split(';')
        for Hit in HitLs:
         if Hit[:4]=='Node':
          Hit=Hit[4:]

          Add='y'
          Tot=0
          Len=len(ExpreLs)
          Ex=0
          HitF='n'
          while Ex<Len:
              Val=ExpreLs[Ex]
              Node=ExpreLs.index[Ex]

              if Node==Hit:

                 HitF='y'
                 if Val==0: Add='n'
              else: Tot+=Val
              Ex+=1
          if Tot!=0: Add='n'

          if Add=='y' and HitF=='y':
 
              In0=list(U.iloc[i])
              In1=In0[:2]+[Hit]+In0[3:]
              In='\t'.join(map(str,In1))+'\n'
              if In not in out:
                 out.append(In)
    GetOut(Out,''.join(out))          
              
                     
def GetOut(OutFile,In):
    OutF=open(OutFile,'w')
    OutF.write(In)
    OutF.close()             
