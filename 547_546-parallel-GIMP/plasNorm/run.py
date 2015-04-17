# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

import os,sys,subprocess

def gnuLive3D(outDict):
   gf=open('gnucomm.xls','w')
   #print>>gf,"set log z"
   print>>gf,"set size square"
   print>>gf,"set xrange[-1:2]"
   print>>gf,"set yrange[-1:2]"
   print>>gf,"set zrange[-1:2]"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   print>>gf,"splot 'history.xls' index 0 using 3:1:2:4 with points pt 5 ps 1 lt palette"
   gf.close()
   subprocess.Popen(['gnuplot','-persist','gnucomm.xls'],shell=False).wait()

def animGif(outDict):
   ni=int(outDict['frameCount'])
   print "Requested Frames:",ni,' ',
   gf=open('gnucomm.xls','w')
   print>>gf,"set term gif animate opt delay 10 size 900,900 x000000"
   print>>gf,"set output 'movie.gif'"
   print>>gf,"set view map"
   print>>gf,"set border linecolor rgb 'cyan'"
   print>>gf,"set cbrange [0:2e10]"
   print>>gf,"set xrange [-.1:1.3]"
   print>>gf,"set yrange [-.2:1.2]"
   print>>gf,"set size square"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,'set title "'+str(i)+' / '+str(ni)+'" tc rgb "cyan"'
      print>>gf,"splot 'history.xls' index "+str(i)+" using 1:2:4 with points pt 5 ps 1.5 lt palette"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def regressionTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'plasNorm-MPM',    'Ncell':4,'shape':'MPM', 'tInt':'momentum','CFL':.5,'load':2e9}]
   ls+=[{'testName':'plasNorm-GIMP',   'Ncell':4,'shape':'GIMP','tInt':'CD','CFL':.5,'load':2e9}]
   ls+=[{'testName':'plasNorm-NewtGM', 'Ncell':4,'shape':'MPM', 'tInt':'dynGM','CFL':4,'load':2e9,'linIter':64,'NewtTol':1e-8,'linTol':1e-6}]
   return ls

if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single
   #dataDict={'mainDir':mainDir,'Ncell':8,'shape':'MPM','tInt':'linCG','CFL':1,'load':2e9,'linTol':1e-3} # not working
   #dataDict={'mainDir':mainDir,'Ncell':8,'shape':'MPM','tInt':'NewtCG','CFL':4,'load':2e9,'NewtTol':1e-5,'linTol':1e-3} # not working
   #dataDict={'mainDir':mainDir,'Ncell':4,'shape':'MPM','tInt':'NewtBCG','CFL':4,'load':2e9,'NewtTol':1e-6,'linTol':1e-3}
   dataDict={'mainDir':mainDir,'Ncell':4,'shape':'MPM','tInt':'dynGM','CFL':4,'load':2e9,'linIter':64,'NewtTol':1e-8,'linTol':1e-6}
   #dataDict={'mainDir':mainDir,'Ncell':4,'shape':'MPM','tInt':'momentum','CFL':.5,'load':2e9}
   #dataDict={'mainDir':mainDir,'Ncell':4,'shape':'GIMP','tInt':'CD','CFL':.5,'load':2e9}
   single(dataDict)
   animGif(dataDict)

