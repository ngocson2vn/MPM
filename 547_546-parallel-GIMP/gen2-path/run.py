# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

import os,sys,subprocess

def tipTrace():
   gf=open('gnucomm.xls','w')
   print>>gf,'set terminal png size 640,480 crop'
   print>>gf,"unset key"
   print>>gf,"set grid"
   #print>>gf,"set xrange[0:1]"
   print>>gf,"set yrange[0:.00375]"
   print>>gf,"set xlabel 'voltage'"
   print>>gf,"set ylabel 'gap'"
   print>>gf,"set output 'probe.png'"
   print>>gf,"plot 'tip.xls' using 1:2 with linespoints pt 5 ps .2"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def gnuLive3D(outDict):
   gf=open('gnucomm.xls','w')
   #print>>gf,"set log z"
   print>>gf,"set size square"
   #print>>gf,"set xrange[-1:2]"
   #print>>gf,"set yrange[-1:2]"
   #print>>gf,"set zrange[-1:2]"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   print>>gf,"splot 'history.xls' index 0 using 1:3:2:4 with points pt 5 ps 1 lt palette",
   print>>gf,",     'vecPlot.xls' index 0 with vectors"
   gf.close()
   subprocess.Popen(['gnuplot','-persist','gnucomm.xls'],shell=False).wait()

def animGif(ni):
   print "Requested Frames:",ni,' ',
   gf=open('gnucomm.xls','w')
   print>>gf,"set term gif animate opt delay 10 size 1600,800 x000000"
   print>>gf,"set output 'movie.gif'"
   print>>gf,"set view map"
   print>>gf,"set border linecolor rgb 'cyan'"
   #print>>gf,"set cbrange [0:100]"
   print>>gf,"set xrange [0:.1]"
   print>>gf,"set yrange [-.004:.008]"
   print>>gf,"set size ratio .3"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,'set title "'+str(i)+' / '+str(ni)+'" tc rgb "cyan"'
      print>>gf,"splot 'nodes.xls'   index "+str(i)+" using 1:2:4 with points pt 5 ps .1 lt 7",
      print>>gf,",     'history.xls' index "+str(i)+" using 1:2:4 with points pt 5 ps 2 lt palette"
      #print>>gf,",     'vecPlot.xls' index "+str(i)+" with vectors lt 5"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def regressionTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'gen2-path-stable','beamL':10,'beamH':1,'beamB':4,'Nthick':1,'maxLoadIter':1,'dtOveride':1./8.,'theta':1,'volt':1500,'linIter':2048,'NewtTol':1e-6,'linTol':1e-8,'np':'4~1~1'}]
   ls+=[{'testName':'gen2-path-pullin','beamL':10,'beamH':1,'beamB':4,'Nthick':1,'maxLoadIter':1,'dtOveride':1./8.,'theta':1,'volt':1600,'linIter':2048,'NewtTol':1e-6,'linTol':1e-8,'np':'4~1~1'}]
   return ls



if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single

   dataDict={'mainDir':mainDir,'Nthick':1,'maxLoadIter':16,'dtOveride':1./32.,'theta':1,'volt':200,'linIter':10000,'NewtTol':1e-6,'linTol':1e-8}

   #dataDict={'mainDir':mainDir,'beamL':.1,'beamB':.05,'Nthick':1,'loadTrend':'linear','maxLoadIter':16,'dtOveride':1./64.,'theta':1,'volt':4000,'linIter':10000,'NewtTol':1e-6,'linTol':1e-8}

   #dataDict={'mainDir':mainDir,'beamL':10,'beamH':1,'beamB':4,'Nthick':1,'dtOveride':1./8.,'theta':1,'volt':1500,'linIter':2048,'NewtTol':1e-6,'linTol':1e-8}
   #dataDict={'mainDir':mainDir,'beamL':10,'beamH':1,'beamB':4,'Nthick':1,'dtOveride':1./8.,'theta':1,'volt':1700,'linIter':2048,'NewtTol':1e-6,'linTol':1e-8}

   single(dataDict)
   tipTrace()
   #gnuLive3D(dataDict)
   #animGif(int(dataDict['frameCount']))
   #animGif(5)











