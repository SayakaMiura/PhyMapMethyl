import os
import numpy
import pandas as pd
import math
from input.Tree import Tree
class Params(object):
    def __init__(self):
        self.message={}
    def parse(self,Argv):
        python='python3'
        Clo2Cell={}
        if Argv[1][-4:]=='.txt':
            Dir=Argv[1].replace(os.sep+Argv[1].split(os.sep)[-1],'')
            Op='Cell' 
            FeaIn0=Argv[1]
        else: 
            Dir=Argv[1] #feature input
            Op='Sample'
            FeaIn0='NA'
        ID=Dir.split(os.sep)[-1]
        print (ID)
        ESLin=''
        Clo2Cell={}
        if Argv[2]=='-Tree':
            NWK=Argv[3]
            if os.path.exists(NWK)!=True: print ('file assigned for -Tree does not exists.',NWK)
            else:
               os.system(python+' MakeIngroup.py '+NWK+' '+Dir)
               ESLin=NWK[:-4]+'_Ingroup.txt'
               print (ESLin)
        if Argv[2]=='-TreeC':
            NWK=Argv[3]
            CellAnno=NWK[:-4]+'.txt'
            if os.path.exists(NWK)!=True: print ('file assigned for -Tree does not exists.',NWK)
            elif os.path.exists(CellAnno)!=True: print ('file assigned for -Tree does not exists.',CellAnno)
            else:
               TrAna=Tree()
               Clo2Cell=TrAna.MakeESLinFromTree(NWK,CellAnno)
               #os.system(python+' MakeIngroup.py '+NWK+' '+Dir)
               ESLin=NWK[:-4]+'_ingroup.txt'
               print (ESLin)
          #  open('a','r').readlines()   
        elif Argv[2]=='-File':
            ESLin=Argv[3]
            if os.path.exists(ESLin)!=True: print ('file assigned for -File does not exists.',ESLin)
        else: print ('tree or list file should be provided')
        ProInfo= Argv[4]
    
      #  if ProInfo=='450': ProInfo='COS450.txt'
      #  elif ProInfo=='RNA': ProInfo='TableS1expressedGene11667COS.txt'

        if os.path.exists(ProInfo)!=True: print ('position file does not exists.',ProInfo)

        print ('Feature ID:',ProInfo)
        return Dir,ID,ESLin,ProInfo,Op,FeaIn0,Clo2Cell
    def GetGroupInfo(self,SysArgv):
        BP2Genels={}
        GeneGroup=''
        GroupLasso='n'
        if SysArgv.count('-Group')!=0 or SysArgv.count('-Group0')!=0:
            GeneGroup='ENTREZ_GENE_ID2UP_KW_BIOLOGICAL_PROCESS.txt'
            if SysArgv.count('-Group')!=0:
              if len(SysArgv)>(SysArgv.index('-Group')+1):
                GeneGroup=SysArgv[SysArgv.index('-Group')+1]
            else: 
              if len(SysArgv)>(SysArgv.index('-Group0')+1):
                GeneGroup=SysArgv[SysArgv.index('-Group0')+1]        
            GeneGroup1=open(GeneGroup,'r').readlines()
            for i in GeneGroup1:
                i=i.split('\t')
                BP=i[1].strip()
                BP2Genels[BP]=BP2Genels.get(BP,[])+[i[0].strip()]   
        if SysArgv.count('-Group')!=0: GroupLasso='y'
        return  BP2Genels,GeneGroup,GroupLasso  
    def GetLambdaSet(self,SysArgv,ESLin):
        Rank='n'
        zLs=[0.1] #change
        yLs=[0]#
        if SysArgv.count('-Group')==0: y=0
        Cut=0
        ID1=ESLin.split(os.sep)[-1][:-4]
        LambdaArg=''
        if SysArgv.count('-Rank')!=0:
            Rank='y'
            Ini1=SysArgv[SysArgv.index('-Rank')+1]
            Fin1=SysArgv[SysArgv.index('-Rank')+2]
            Ini2=SysArgv[SysArgv.index('-Rank')+3]
            Fin2=SysArgv[SysArgv.index('-Rank')+4]
            LambdaArg='--initial_lambda1 '+Ini1+' --final_lambda1 '+Fin1+' --initial_lambda2 '+Ini2+' --final_lambda2 '+Fin2
            ID1=ESLin.split(os.sep)[-1][:-4]
        if SysArgv.count('-Range')!=0:
            Ini1=SysArgv[SysArgv.index('-Range')+1]
            Fin1=SysArgv[SysArgv.index('-Range')+2]
            Ini2=SysArgv[SysArgv.index('-Range')+3]
            Fin2=SysArgv[SysArgv.index('-Range')+4].strip()
            if Ini1==Fin1: zLs=[float(Ini1)]
            else: 
               # zLs= list(numpy.arange(float(Ini1),float(Fin1),0.1))
                zLs= numpy.logspace(math.log10(float(Ini1)),math.log10(float(Fin1)), 3, endpoint=True)
            if Ini2==Fin2: yLs=[float(Ini2)]
            else: 
              #  yLs=list(numpy.arange(float(Ini2),float(Fin2),0.1))         
                yLs= numpy.logspace(math.log10(float(Ini2)),math.log10(float(Fin2)), 3, endpoint=True)                    
        return LambdaArg,Rank,zLs,yLs,ID1  
    def GetGeneinInfo(self,ProInfo):
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
    def GetTuOrder(self,ESLRes):
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
    def MakeMapIn(self,ProbNum,ProbLs1,MapIn):     
       Map=['0\t\n']
       c=0
       while c<ProbNum:
           Map.append(str(c+1)+'\t'+ProbLs1[c].strip()+'\n')
           c+=1
       self.GetOut(MapIn,''.join(Map).replace(' ',''))         
    def MakeFeatureFilePreprocess(self,Dir,TuOr):
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
    def MakeESLInputFiles_simple(self,ESLin,Dir,ID1):
     # print (ESLin)
      ESLRes=pd.read_csv(ESLin,sep='\t') #tree.txt

      TuOrIn=Dir+os.sep+ID1+'_feature_input_sample.txt'
      FeaIn=Dir+os.sep+ID1+'_feature_input.txt'
      MapIn=Dir+os.sep+ID1+'_feature_mapping_input.txt' 
      GroIn=Dir+os.sep+ID1+'_group_indices_input.txt'    
      FeildIn=Dir+os.sep+ID1+'_field_input.txt'
      if os.path.exists(FeaIn)==True:
          TuOr0=open(TuOrIn,'r').readlines()
          TuOr=[]
          for i in TuOr0:
              TuOr.append(i.strip())
      else:
        TuOr=self.GetTuOrder(ESLRes)
        self.GetOut(TuOrIn,'\n'.join(TuOr))
        CellOr=[]
        Amat,FileLs,ProbLs=self.MakeFeatureFilePreprocess(Dir,TuOr)

        PathIn=Dir+os.sep+'path.txt'
        self.GetOut(PathIn,'.fas\n'.join(ProbLs)+'.fas\n') 
    
        if len(TuOr)==len(FileLs):
           print ('remove probes with nan. Making table of beta values...')
           FeaTa=pd.DataFrame(Amat,columns=ProbLs)#pd.read_csv(FeaIn0,sep=',',names=ProbLs)
           print ('pruning')
           FeaTaClean=FeaTa.dropna(axis=1)
           print (FeaTaClean)
           A=FeaTaClean.max()
           if max(list(A))>1:
               print ('features are >1 so normalize it')
               FeaTaClean = (FeaTaClean-FeaTaClean.min())/(FeaTaClean.max()-FeaTaClean.min())
           print (FeaTaClean)
          # open('a','r').readlines()
           FeaTaClean.to_csv(FeaIn,sep=',',header=False,index=False)
           ProbLs1=list(FeaTaClean.columns.values)
           Go='y'
        else: Go='n'   
           
        if Go=='y':    
           ProbNum=len(ProbLs1)   
           print (ProbNum,len(ProbLs1),ProbLs1[0],ProbLs1[-1])
         #  open('a','r').readlines()
    #   MapIn=Dir+os.sep+ID1+'_feature_mapping_input.txt' 
           self.MakeMapIn(ProbNum,ProbLs1,MapIn)     
           self.GetOut(Dir+os.sep+ID1+'_feature_input_gene.txt','\n'.join(ProbLs1)+'\n')
           PosLs=list(range(1,ProbNum+1))
           NumStr=','.join(map(str,PosLs))
           Seq=list(range(1,ProbNum+1))
           Seq1=list(map(str,Seq))
           PosAllSep=','.join(Seq1)
           AllOne=['1']*ProbNum
           AllOne1=','.join(AllOne)
           print (len(AllOne),ProbNum,len(Seq),len(Seq1),PosAllSep.count(','),AllOne1.count(','))

           #    print (GroIn)
         #  GetOut(GroIn,PosAllSep+'\n'+PosAllSep+'\n'+AllOne1+'\n')
           self.GetOut(GroIn,'1\n'+str(ProbNum)+'\n'+str(ProbNum**0.5)+'\n')  
      # open('a','r').readlines() 
      

      return ESLRes,TuOr,FeaIn,GroIn,MapIn, FeildIn                       
    def MakeESLInputFiles(self,ESLin,Dir,BP2Genels,GroupLasso,Entr2Pro,ID1,Op,FeaIn0,ID):
      print (ESLin)
      ESLRes=pd.read_csv(ESLin,sep='\t')
     # print (ESLRes)
     # open('a','r').readlines()
      if Op=='Cell':
         TuOr=self.GetTuOrder(ESLRes)
      else: 
          if ID[0]=='R': TuOrIn=Dir+os.sep+ID+'_feature_input_sample.txt'
          else: TuOrIn=Dir+os.sep+'input_feature_input_sample.txt'
          TuOr0=open(TuOrIn,'r').readlines()
          TuOr=[]
          for i in TuOr0:
              TuOr.append(i.strip())         
      FeaIn=Dir+os.sep+ID1+'_feature_input.txt'
      MapIn=Dir+os.sep+ID1+'_feature_mapping_input.txt' 
      GroIn=Dir+os.sep+ID1+'_group_indices_input.txt'    
      FeildIn=Dir+os.sep+ID1+'_field_input.txt'
      if os.path.exists(FeaIn)!=True:
    #   if Op!='Cell' and Op!='Sample': 
    #      print ('Op is Sample or Cell',Op)
    #      open('a','r').readlines()
    #   else:#elif Op=='Cell': 
          # print (FeaIn0)
           FeaIn=FeaIn0

           ProbLs0=open(FeaIn0[:-4]+'_gene.txt','r').readlines()
           ProbLs1=[]
           for i in ProbLs0: ProbLs1.append(i.strip())
          # CellOr=[]
      #     Go='y' 
           
       #if Go=='y':    
           ProbNum=len(ProbLs1)   
           print (ProbNum,len(ProbLs1),ProbLs1[0])
    #   MapIn=Dir+os.sep+ID1+'_feature_mapping_input.txt' 
           self.MakeMapIn(ProbNum,ProbLs1,MapIn)     
       
           PosLs=list(range(1,ProbNum+1))
           NumStr=','.join(map(str,PosLs))
           Seq=list(range(1,ProbNum+1))
           Seq1=list(map(str,Seq))
           PosAllSep=','.join(Seq1)
           AllOne=['1']*ProbNum
           AllOne1=','.join(AllOne)
           print (len(AllOne),ProbNum,len(Seq),len(Seq1),PosAllSep.count(','),AllOne1.count(','))
           if BP2Genels!={}:
               BP2genels=GetBP2genels(BP2Genels,Entr2Pro,ProbLs)
                     
           if GroupLasso=='n':#BP2Genels=={}: 
               print (GroIn)
         #  GetOut(GroIn,PosAllSep+'\n'+PosAllSep+'\n'+AllOne1+'\n')
               self.GetOut(GroIn,'1\n'+str(ProbNum)+'\n'+str(ProbNum**0.5)+'\n')  
      # open('a','r').readlines() 
      
           else:  
               Gdf={}
               ProPosOrAll=[]
               Minor=[]
               BP2genels={}
               AllBPgeneLs=[]
               S=1
               for BP in BP2Genels:
                   Genels=BP2Genels[BP]
                   Genels1=[]
                   ProPosOr=[]
           
                   for G in Genels:
                       Pro=Entr2Pro.get(G,'')
               
                       if Pro in ProbLs:
               #    print (Pro,ProbLs.index(Pro))
                           ProPosOr.append(ProbLs.index(Pro)+1)
                           BP2genels[BP]=BP2genels.get(BP,[])+[Pro]
                           AllBPgeneLs.append(Pro)
                   if len(ProPosOr)<2: 
                      Minor+=ProPosOr
                      print ('<2',BP,Genels)
              #open('a','r').readlines()
                   else: 
                       ProPosOrAll+=ProPosOr
                       E=S+len(ProPosOr)-1
                       Gdf[BP]=[str(S),str(E),(float(len(ProPosOr)))**0.5]
                       S=E+1    
                #   open('a','r').readlines()
               Miss=[]
               c=1
               while c<=len(ProbLs):
                   if c not in ProPosOrAll: Miss.append(c)
                   c+=1
               print (len(Miss),c)
       #open('a','r').readlines()    
               ProPosOrAll+=Miss    
               E=S+len(Miss)-1    
               Gdf['Miss']=[str(S),str(E),(float(len(Miss)))**0.5] 
               GetOut(FeildIn,','.join(map(str,map(int,ProPosOrAll))))
               Gdf1=pd.DataFrame(Gdf)
               print (Gdf1)
               Gdf1.to_csv(GroIn[:-4]+'_withHeader.csv',index=False)
               Gdf1.to_csv(GroIn,index=False,header=False)    
      return ESLRes,TuOr,FeaIn,GroIn,MapIn, FeildIn
    def InvertDic2(self,St2Seq):
        Hap2ID={}
        for St in St2Seq:
          SeqLs=St2Seq[St]
          for Seq in SeqLs:
            Hap2ID[Seq]=St
        return Hap2ID        
    def GetTuCellLs(self,Clo2Cell,CellLsFile,Op):
        CellOr=[]
        Cell2Tu={}   
        if Op=='Cell':
             CellLs=open(CellLsFile,'r').readlines()
             for i in CellLs: 
                 CellOr.append(i.strip())
             Cell2Tu=self.InvertDic2(Clo2Cell)             
         #    CellClo=open(CellCloTaFile,'r').readlines()[1:]
          #   for i in CellClo:
           #      i=i.split('\t')
            #     Cell=i[0].strip()
             #    if Cell in CellOr:
              #       Cell2Tu[Cell]=i[1].strip()
  
        return CellOr,Cell2Tu  
    def MakeResponseFile(self,CellLs,Cell2Tu,TuOr,Ingroup,Outgroup,Op,Dir,IDn,ID1):
    # print (CellLs[0])
     Res=[]
     Res1=[]
     PosC=0
     NegC=0
     if Op=='Sample':
       for Tu in TuOr:
          if Tu in Ingroup:
              Res.append('1\n')
              Res1.append(Tu+'\t1\n')
              PosC+=1
          elif Tu in Outgroup: 
              Res.append('-1\n')
              Res1.append(Tu+'\t-1\n')
              NegC+=1
          #else:
          #   print (Tu) 
          #   open('a','r').readlines()
     elif Op=='Cell':
        print ('making response file')
       # print (Ingroup,Outgroup)
        for Cell in CellLs:
            Tu=Cell2Tu[Cell]
            if Tu in Ingroup:
              Res.append('1\n')
              Res1.append(Cell+'\t1\n')
              PosC+=1
            elif Tu in Outgroup:
              Res.append('-1\n')
              Res1.append(Cell+'\t-1\n')
              NegC+=1
            elif Tu=='Outgroup': print ('Outgroup cell annotation',Cell)  
            else:
              print (Cell,Tu) 
              open('a','r').readlines()     
     else:                                   
        print ('Op should be Cell or Sample',Op) 
        open('a','r').readlines()        
     ResIn= Dir+os.sep+IDn+'_'+ID1+'_response_input.txt'   
     ResIn1=Dir+os.sep+IDn+'_'+ID1+'_response_input1.txt'
     SweiIn=Dir+os.sep+IDn+'_'+ID1+'_sweights_input.txt'
     self.GetOut(ResIn,''.join(Res))
     self.GetOut(ResIn1,''.join(Res1))
   #  print (PosC,NegC,Outgroup,Ingroup)#,Cell2Tu)
     if NegC>0 and PosC>0:
         PosNeg=1.0*PosC/NegC
         self.GetOut(SweiIn,str(PosNeg)+'\n1.0')  
     return ResIn,ResIn1,SweiIn  
    def  GetInfo(self,ResIn,TuOr):                 
         ResOr=open(ResIn,'r').readlines()
         SampHead=[]
         c=0
         while c<len(ResOr):
                     SampHead.append(ResOr[c].strip()+'_'+TuOr[c])
                     c+=1
       #  print (SampHead) 
       
         return ResOr,SampHead   
    def GetOut(self,OutFile,out):
        OutF=open(OutFile,'w')
        OutF.write(out)
        OutF.close()                         
