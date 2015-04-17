# Philip Wallstedt 2004-2010

import os,sys,subprocess

def gnuLive3D(outDict):
   gf=open('gnucomm.xls','w')
   #print>>gf,"set log z"
   #print>>gf,"set size square"
   #print>>gf,"set xrange[-1:2]"
   #print>>gf,"set yrange[-1:2]"
   #print>>gf,"set zrange[-1:2]"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   print>>gf,"splot 'history.xls' index 0 using 3:1:2:4 with points pt 5 ps 1 lt palette"
   gf.close()
   subprocess.Popen(['gnuplot','-persist','gnucomm.xls'],shell=False).wait()

def regressionTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'massMPI','Ncell':8,'shape':'MPM','tInt':'momentum','CFL':.9,'np':'2~2~2'}]
   return ls



if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single
   #dataDict={'mainDir':mainDir,'Ncell':8,'shape':'MPM','tInt':'mass','CFL':.9,'np':'1~1~1'}
   #dataDict={'mainDir':mainDir,'Ncell':4,'shape':'MPM','tInt':'mass','CFL':.9,'np':'2~1~1'}
   dataDict={'mainDir':mainDir,'Ncell':8,'shape':'GIMP','tInt':'CD','CFL':.9,'np':'2~2~2'}
   single(dataDict,False)
   #for el in outDict:print el,'=',outDict[el]
   gnuLive3D(dataDict)



