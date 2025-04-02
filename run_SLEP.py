import os
import sys
import subprocess


Out=sys.argv[1]
WeiTa=Out+'.txt'
FeaIn=sys.argv[2]
z=Out.split('_feature_weight')[-1].split('_')[0]
y=Out.split('_feature_weight')[-1].split('_')[1]
SweiIn=sys.argv[3]
GroIn=sys.argv[4]
ResIn=sys.argv[5] 
def GetOut(Out,out):
    OutF=open(Out,'w')
    OutF.write(out)
    OutF.close()
Com='./sg_lasso -f '+FeaIn+' -z '+str(z)+' -y '+str(y)+' -s '+SweiIn+' -n '+GroIn+' -r '+ResIn+' -w '+Out
if os.path.exists(WeiTa)!=True:
  print ('doing SLEP...')
  Cwd=os.getcwd()
  os.chdir('bin')
  os.system(Com)
  os.chdir(Cwd)  
  if os.path.exists(Out+'.xml')==True:
   print ('finished SLEP. begin grep')

   Com='grep -P \"<item>.*</item>\" '+Out+'.xml | sed -re \"s/.*<item>(.*)<\\/item>.*/\\1/\" > '+Out+'_temp.txt'
   os.system(Com)
   print ('finished grep. begin paste')

   MapIn=sys.argv[6]
   Com='paste <(sed -e \"1d\" '+MapIn+') '+Out+'_temp.txt | grep -v \"0.00000000000000000e+00\" > '+WeiTa+'\n'
   if os.path.exists('Test.sh')==True: os.remove('Test.sh')
   GetOut('Test.sh','#!/bin/bash\n'+Com)
   os.system('chmod +x '+os.getcwd()+os.sep+'Test.sh')
   subprocess.call(os.getcwd()+os.sep+'Test.sh')
   print ('finished paste')

