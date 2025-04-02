import os
import pandas as pd
import glob
import sys
#HitFeaNodeTa='E:\\Desktop\\RNA\\Test\\C537_a_Tree_lambda10.1_lambda20.1_out_feature_weight_Hit.txt'
Nwk=sys.argv[1]#'E:\\Desktop\\Methylation\\HCC8010\\TreeMethyl1NoN.nwk'
HitFeaNodeTa=Nwk[:-4]+'_HitGeneNode.txt'#'E:\\Desktop\\Methylation\\HCC8010\\TreeMethyl1NoN_HitGeneNode0.2.txt'
#Nwk='E:\\Desktop\\RNA\\Test\\C537edit.nwk'

NodeIDin=Nwk[:-4]+'.txt'
#os.system('Rscript --vanilla Runl1ouWOboo.r '+HitFeaNodeTa+' '+Nwk) ###need
#open('a','r').readlines()
print ('summarize')
#Boo=HitFeaNodeTa+'_boo.txt'
Hit=HitFeaNodeTa+'_Shift.txt'
NodeID=HitFeaNodeTa+'_TreeNode.txt'
NodeID=pd.read_csv(NodeID,sep='\t')
Hit=pd.read_csv(Hit,sep='\t')
In=pd.read_csv(HitFeaNodeTa,sep='\t')
#Boo=pd.read_csv(Boo,sep=',')
NodeIDin=open(NodeIDin,'r').readlines()[1:]
NodeID2TipLsIn={}
for i in NodeIDin:
    i=i.split('\t')
   # print (i)
    ID=i[0]
    Ls=i[1].split(';')
    In0=[]
    for ii in Ls:
        ii=ii[(len(ii.split('-')[0])+1):]
        In0.append(ii)
    In0.sort()    
    NodeID2TipLsIn[ID]=In0 
print (NodeID2TipLsIn)
#open('a','r').readlines()    


def GetNodeID2TipLs(Ta,Gene):
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
            TipLs.sort()
            Node2TipLs[N]=TipLs
   # print (Node2TipLs)
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
c=0

RowN=len(In)
outS='y'
OutUnS='y'
print (In)
while c<RowN:
    Gene=In["Gene"].tolist()[c]
    NodeIDInLs=In["Node"].tolist()[c].split(';')
    NodeID2TipLs=GetNodeID2TipLs(NodeID,Gene)
    HitBra=Hit[Hit['Gene']==Gene]
 #   BooVal=Boo[Boo['Unnamed: 0']==Gene]
 #   print (BooVal)
   # print(NodeID2TipLs)
    if len(HitBra)>0:
        HitBraLs=HitBra['branch'].tolist()
      #  print (HitBraLs)
        Entered='n'
        for Bra in HitBraLs:
            TipLs=NodeID2TipLs['X'+str(Bra)]
            for NodeIDin in NodeIDInLs:
                TipLsIn=NodeID2TipLsIn[NodeIDin]
              #  print (Bra,TipLs)
                Ori0=In[In["Gene"]==Gene].to_dict()
                Add=AddColumn(Ori0,{'Node':[NodeIDin]}) #dataframe              
                if TipLs==TipLsIn:
                    if outS=='y':
                        if Entered=='n': 
                            out=In[In["Gene"]==Gene]
                        Entered='y'
                        outS='n'  
                        out1=Add
                    else:
                        if Entered=='n': out=pd.concat([out,In[In["Gene"]==Gene]])
                        Entered='y'
                        out1=pd.concat([out1,Add])
                  #  print (out)  
                else:

              
                    if OutUnS=='y':

                        outUn=Add#In[In["Gene"]==Gene]
                        OutUnS='n'

                    else:
                    
                        outUn=pd.concat([outUn, Add])
                      #  print (outUn)
                     #   open('a','r').readlines()                        
                   # print (out)                  
      #  open('a','r').readlines()
    c+=1
outUn.to_csv(HitFeaNodeTa[:-4]+'_UnSupported.txt',sep='\t', index=False) 
out1.to_csv(HitFeaNodeTa[:-4]+'_Supported.txt',sep='\t', index=False)   
HitTaFil=HitFeaNodeTa[:-4]+'_Fil.txt'      
Boo='n'
if Boo=='y':  
 
 out.to_csv(HitTaFil,sep='\t', index=False)    
 os.system('Rscript --vanilla Runl1ou.r '+HitTaFil+' '+Nwk)   ##need
#open('a','r').readlines()
 print ('summarize bootstrap res')
 HitFeaNodeTa=HitTaFil
 Boo0=HitFeaNodeTa+'_boo.txt'
 tmp=open(Boo0,'r').readlines()
 H=['Gene,'+tmp[0]]+tmp[1:]
 Boo=Boo0[:-4]+'1.txt'
 OutF=open(Boo,'w')
 OutF.write(''.join(H))
 OutF.close()
 Hit=HitFeaNodeTa+'_Shift.txt'
 NodeID=HitFeaNodeTa+'_TreeNode.txt'
 NodeID=pd.read_csv(NodeID,sep='\t')
 Hit=pd.read_csv(Hit,sep='\t')
 #In=pd.read_csv(HitFeaNodeTa,sep='\t')
 Boo=pd.read_csv(Boo,sep=',') 
#print (Boo)

 c=0
 out=''
 S='y'
 while c<RowN:
    Gene=In["Gene"].tolist()[c]
    NodeIDInLs=In["Node"].tolist()[c].split(';')
    NodeID2TipLs=GetNodeID2TipLs(NodeID,Gene)
    HitBra=Hit[Hit['Gene']==Gene]
    BooVal=Boo[Boo['Gene']==Gene]
  #  print (BooVal)
 #   print(NodeID2TipLs)
    if len(HitBra)>0:
        HitBraLs=HitBra['branch'].tolist()
   #     print (HitBraLs)
        
        for Bra in HitBraLs:
            TipLs=NodeID2TipLs['X'+str(Bra)]
            for NodeIDin in NodeIDInLs:
                TipLsIn=NodeID2TipLsIn[NodeIDin]
               # print (Bra,TipLs)
                if TipLs==TipLsIn:
                    BooValTar=BooVal['V'+str(Bra)]
                   # print (BooValTar)
                    Ori0=In[In["Gene"]==Gene]
                    Add=AddColumn(Ori0,{'BooSup':[BooValTar[list(BooValTar.keys())[0]]],'Node':[NodeIDin]}) #dataframe
                    if S=='y':
                        out=Add
                        S='n'
                    else:
                        out=pd.concat([out, Add])
                  #  print (out)                        
       # open('a','r').readlines()
    c+=1

 out.to_csv(HitFeaNodeTa[:-8]+'_Supported.txt',sep='\t', index=False)  
os.system('Rscript --vanilla Heatmap.r '+HitFeaNodeTa[:-4]+'_Supported.txt')
RmLs=glob.glob(HitFeaNodeTa[:-8]+'.txt_*')+glob.glob(HitFeaNodeTa[:-8]+'_Fil*')
for i in RmLs:
   os.remove(i)