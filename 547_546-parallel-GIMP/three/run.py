# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

import os,sys,subprocess

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
   ls+=[{'testName':'three-momGM','Nmul':1,'tInt':'momGM','CFL':128,'load':.001,'mode':1,'linIter':1024,'NewtTol':1e-6,'linTol':1e-8}]
   return ls


if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single

   dataDict={'mainDir':mainDir,'Nmul':1,'tInt':'momGM','CFL':128,'load':.001,'mode':1,'linIter':1024,'NewtTol':1e-6,'linTol':1e-8}
   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'dynGM','CFL':128,'load':1e-5,'mode':1,'linIter':256,'NewtTol':1e-6,'linTol':1e-8}
   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'momentum','CFL':.9,'load':1e-5,'mode':1}

   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'NewtCG','CFL':64,'load':1e-5,'mode':1,'linIter':400,'linTol':1e-8}


   single(dataDict,False)
   #print 'wall time',dataDict['wallTime']
   tipTrace(dataDict)











