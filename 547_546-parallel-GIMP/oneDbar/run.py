# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

import os,sys,subprocess

def probeTrace(outDict):
   gf=open('gnucomm.xls','w')
   print>>gf,'set terminal png size 640,480 crop'
   print>>gf,"set key bottom left"
   print>>gf,"set grid"
   print>>gf,"set xlabel 'Time'"
   #print>>gf,"set ylabel 'Displacement'"
   print>>gf,"set ylabel 'Stress'"
   print>>gf,"set output 'probeTrace.png'"
   print>>gf,"plot 'history.xls' using 1:5 title '    L' with linespoints pt 5 ps .2",
   print>>gf,",    'history.xls' using 1:4 title '3/4 L' with linespoints pt 5 ps .2",
   print>>gf,",    'history.xls' using 1:3 title '1/2 L' with linespoints pt 5 ps .2",
   print>>gf,",    'history.xls' using 1:2 title '1/4 L' with linespoints pt 5 ps .2"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()



if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single

   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':2,'Nwide':1,'shape':'MPM','tInt':'NewtCG','CFL':128,'load':1e-5,'mode':1,'linIter':400,'NewtTol':1e-6,'linTol':1e-8}
   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'dynGM','CFL':128,'load':1e-5,'mode':1,'linIter':256,'NewtTol':1e-6,'linTol':1e-8}
   dataDict={'mainDir':mainDir,'Nlong':100,'CFL':.9,'load':1e-5}
  
   single(dataDict)
   probeTrace(dataDict)











