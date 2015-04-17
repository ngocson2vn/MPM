# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

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
   print>>gf,"splot 'history.xls' index 0 using 1:3:2:4 with points pt 5 ps 1 lt palette"
   gf.close()
   subprocess.Popen(['gnuplot','-persist','gnucomm.xls'],shell=False).wait()

def tipTrace(outDict):
   gf=open('gnucomm.xls','w')
   print>>gf,'set terminal png size 640,480 crop'
   print>>gf,"set key top left"
   print>>gf,"set grid"
   print>>gf,"set xlabel 'Time'"
   print>>gf,"set ylabel 'Velocity'"
   print>>gf,"set output 'tipTrace.png'"
   print>>gf,"plot 'tip.xls' using 1:2 with linespoints pt 5 ps .2",
   print>>gf,",    'tip.xls' using 1:3 with linespoints pt 5 ps .2"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def regressionTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'fixed-uniform-quasiCG','beamL':.1,'Nthick':1,'Nwide':1,'tInt':'quasiCG','dtOveride':.3,'load':1,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8,'np':'4~1~1'}]
   ls+=[{'testName':'fixed-uniform-quasiGM','beamL':.1,'Nthick':1,'Nwide':1,'tInt':'quasi'  ,'dtOveride':.3,'load':1,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8,'np':'3~1~1'}]
   return ls


if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single
   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'tInt':'quasiCG','dtOveride':.3,'load':1,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8}
   dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'tInt':'quasiCG','dtOveride':.3,'load':1,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8,'np':'2~1~1'}
   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'quasi'  ,'dtOveride':.3,'load':1,'linIter'  :1000,'NewtTol':1e-6,'linTol':1e-8}
   single(dataDict,False)
   tipTrace(dataDict)
   #gnuLive3D(dataDict)











