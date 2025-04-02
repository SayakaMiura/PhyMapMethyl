import os
import glob
import sys
import shutil
import pandas as pd

#python3 PhyMapMethyl.py /home/esl/Desktop/PhyMap/Example/input /home/esl/Desktop/PhyMap/Example/Tree.nwk Tree5 
#python3 PhyMapMethyl.py [path to input beta value files] [path to tree file] "ID" 

Fea=sys.argv[1]
Tree=sys.argv[2]
ID=sys.argv[3]
Gene='NA'


if os.path.exists(Tree)!=True:
    print (Tree, "file does not exists. Try it again.")
else:
    InputType='-Tree'
    WeiCut=[0.001,0.0001]
    MaxSigleHit=10000#1000
    MinSigleHit=1000#100
    RepMax=[5,30,40,50]#[1,3,5,7,10]
    IteMax=-1#int(sys.argv[4]) #2
    FPcut=0.1#float(sys.argv[5]) #0.5
    FPcutMax=0.1  
    Py='python3 RunESL10.py'  

    print ('Making clades file...')
    os.system('python MakeIngroup1.py '+Tree+' '+Fea+' '+ID) 
 
def Del(Ls):
   for i in Ls:
       if os.path.exists(i)==True: os.remove(i)  
def GetFP(Adf,Sdf):
    A=Adf['Gene'].to_list()
    A=len(list(set(A)))
    S=Sdf['Gene'].to_list()
    S=len(list(set(S)))
    if A==0: FP=0
    else:
       FP=1.0*(A-S)/A
    print ('all hit and supported hit: ',A,S,FP)
    return A,S,FP  
def Check(Note,IteMax,I):
    iup=I  
    FPnew=0
    for i in Note:
        if i.find('terminated')!=-1: iup= IteMax+1
        elif i.find('There is no shift')!=-1: iup= IteMax+1  
        if i.find('FP rate:')!=-1: FPnew=float(i.split(':')[-1])
    return iup,FPnew                
if Fea[-4:]!='.txt': Dir=Fea+'/'
else: Dir='/'.join(Fea.split(os.sep)[:-1])+'/'
print ('output directory: ',Dir)       
Ite='y'
W=0
R=0
I=0
Log=[]
FPnewPre=0
while Ite=='y':
    print ('iteration ',I, ' Prepare the analysis...')

    PreFileLs=glob.glob(Dir+'*_AllHit*.txt')
    if I==0:
       PreFileLs+=glob.glob(Dir+'PhyMap_*.txt')+glob.glob(Dir+'*_HitGeneNode.txt')
    if R==0:
        PreFileLs+=glob.glob(Dir+'*_HitGeneNode.txt')
        directory_path = Dir+'Rep1'
        if os.path.exists(directory_path)==True:
             shutil.rmtree(directory_path)
             
    print ('Delete output files in the output directory',len(PreFileLs))  
    Del(PreFileLs)                  
    Arg=' '.join(map(str,[Py,Fea,InputType,Tree,WeiCut[W],RepMax[R],MaxSigleHit,MinSigleHit,Gene]))
    print (Arg)
    os.system(Arg)
    print (Arg)
    print ('iteration ',I, ' end the analysis.')
    Sfile=Dir+Tree.split(os.sep)[-1][:-4]+'_AllHitNode1_Supported.txt'
    if os.path.exists(Sfile)!=True:
        Sfile=Dir+Tree.split(os.sep)[-1][:-4]+'_AllHitNode11_Supported2.txt'
    Afile=Dir+Tree.split(os.sep)[-1][:-4]+'_AllHitNode.txt'
    if os.path.exists(Afile)==True and os.path.exists(Sfile)==True:
     if I==0:
          Sdf=pd.read_csv(Sfile,sep='\t')
          Adf=pd.read_csv(Afile,sep='\t')
          Sdf.to_csv(Dir+'PhyMap_Supported.txt',sep='\t', index=False)
          Adf.to_csv(Dir+'PhyMap_All.txt',sep='\t', index=False) 
     else:
          Adds=pd.read_csv(Sfile,sep='\t')
          Adda=pd.read_csv(Afile,sep='\t')    
          Sdf=pd.concat([Sdf, Adds]) 
          Adf=pd.concat([Adf, Adda])  
          Sdf.to_csv(Dir+'PhyMap_Supported.txt',sep='\t', index=False)
          Adf.to_csv(Dir+'PhyMap_All.txt',sep='\t', index=False)   
           
     A,S,FP=GetFP(Adf,Sdf)

     Note=glob.glob(Dir+Tree.split(os.sep)[-1][:-4]+'*_Note.txt')
     Note=open(Note[0],'r').readlines()           

     Log.append('\nOverall all hit, supported, and FP: '+str(A)+' '+str(S)+' '+str(FP)+'\n')              
    else: 
      
      Note=['SLEP l1ou pipeline did not detect shifts.']                 
    Log+=['Iteration '+str(I)+'\n'+Arg+'\n']+Note
    I,FPnew=Check(Note,IteMax,I)
    if I>IteMax or FPnew>FPcutMax: Ite='n'
    else:
       DeltaFP=FPnew-FPnewPre
       if DeltaFP<0.1 and FPnew>FPcut: Ite='n'
       else:
           print ('Too little FP. So, repeat...')

           if I>IteMax: Ite='n'
           W+=1
           R+=1
           I+=1    
    print (''.join(Log))
    FPnewPre=FPnew 

OutF=open(Dir+'log.txt','w')
OutF.write(''.join(Log))
OutF.close()   
Sup=Dir+'PhyMap_Supported.txt'  
if InputType=='-Tree':
        Com='Rscript --vanilla Heatmap3.r '+Sup
        os.system(Com)
        print (Com)

      
os.remove(Dir+'path.txt')
    
DelLs=glob.glob(Dir+'*_out_feature_weight*')  
DelLs+=glob.glob(Dir+'input_*.txt') 
print (Tree.split(os.sep)[-1][:-4]) 
DelLs+=glob.glob(Dir+Tree.split(os.sep)[-1][:-4]+'_*.txt')  

for i in DelLs:
   if os.path.exists(i)==True:
    os.remove(i)
directory_path = Dir+'Rep1'
if os.path.exists(directory_path)==True:
     shutil.rmtree(directory_path)  
os.remove(Dir+'log.txt')
os.remove(Dir+'PhyMap_All.txt')
