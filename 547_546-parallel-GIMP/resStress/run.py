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
   #print>>gf,"set yrange [-.0001:.0001]"
   print>>gf,"set xlabel 'load factor'"
   print>>gf,"set output 'deflection.png'"
   print>>gf,"set ylabel 'deflection'"
   print>>gf,"plot 'tip.xls' using 1:2 title 'mid' with linespoints pt 5 ps .5",
   print>>gf,",    'tip.xls' using 1:3 title 'edge' with linespoints pt 4 ps 1"
   #print>>gf,",    'tip.xls' using 1:4 title 'exact' with linespoints pt 3 ps 1"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def exactForce():
   E=158.1e3
   L=.4
   b=.12
   h=.004
   d=.0001
   P=16.*E*b*h*h*h*d/(L*L*L) # (N)
   area=.002*.03
   print 'exact force',P
   print 'exact pressure',P/area

def regressionTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'resStress-good','beamL':.1,'beamB':.024,'Nthick':2,'Nwide':4,'dtOveride':.2,'fM':2e-3,'RS':400,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8,'np':'4~1~1'}]
   ls+=[{'testName':'resStress-fail','beamL':.1,'beamB':.022,'Nthick':2,'Nwide':4,'dtOveride':.2,'fM':2e-3,'RS':400,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8,'np':'4~1~1'}]
   return ls


if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single
   dataDict={'mainDir':mainDir,'beamL':.1,'beamB':.024,'Nthick':2,'Nwide':4,'dtOveride':.4,'fM':2e-3,'RS':0,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8} # succeeds
   #dataDict={'mainDir':mainDir,'beamL':.1,'beamB':.024,'Nthick':2,'Nwide':4,'dtOveride':.2,'fM':2e-3,'RS':400,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8} # succeeds
   #dataDict={'mainDir':mainDir,'beamL':.1,'beamB':.022,'Nthick':2,'Nwide':4,'dtOveride':.2,'fM':2e-3,'RS':400,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8} # fails
   #dataDict={'mainDir':mainDir,'beamL':.1,'beamB':.024,'Nthick':2,          'dtOveride':.2,'fM':2e-3,'RS':400,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8} # succeeds
   #dataDict={'mainDir':mainDir,'beamL':.1,'beamB':.022,'Nthick':2,          'dtOveride':.2,'fM':2e-3,'RS':400,'linIter':1000,'NewtTol':1e-6,'linTol':1e-8} # fails
   #exactForce()
   single(dataDict)
   tipTrace(dataDict)
   #gnuLive3D(dataDict)











