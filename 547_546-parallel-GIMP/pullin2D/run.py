# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

import os,sys,subprocess

def tipTrace():
   gf=open('gnucomm.xls','w')
   print>>gf,'set terminal png size 640,480 crop'
   print>>gf,"set key bottom left"
   print>>gf,"set grid"
   print>>gf,"set yrange [0:.004]"
   print>>gf,"set xlabel 'Volt'"
   print>>gf,"set ylabel 'Gap'"
   print>>gf,"set output 'pullin.png'"
   print>>gf,"plot 'tip.xls'      using 1:2 title 'MPM' with linespoints lw 2 pt 5 ps .5",
   print>>gf,",    'pullinEB.xls' using 1:2 title 'EB'  with lines lw 2"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def animGif(ni):
   ni=int(ni)
   print "Requested Frames:",ni,' ',
   gf=open('gnucomm.xls','w')
   print>>gf,"set term gif animate opt delay 100 size 1600,400 x000000"
   print>>gf,"set output 'movie.gif'"
   print>>gf,"set view map"
   print>>gf,"set border linecolor rgb 'cyan'"
   #print>>gf,"set cbrange [0:100]"
   print>>gf,"set xrange [0:.4]"
   print>>gf,"set yrange [0:.004]"
   print>>gf,"set size ratio .25"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,'set title "'+str(i)+' / '+str(ni)+'" tc rgb "cyan"'
      print>>gf,"splot 'history.xls' index "+str(i)+" using 1:2:4 with points pt 5 ps 1 lt palette"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def regressionTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'pullin2D-geom-3000','beamL':.1,'pois':.31,'Nthick':1,'lastFac':.01,'Nstep':8,'theta':1,'volt':3000,'linIter':200,'NewtTol':1e-6,'linTol':1e-8,'maxLoadIter':32,'np':'4~1~1'}]
   ls+=[{'testName':'pullin2D-geom-3020','beamL':.1,'pois':.31,'Nthick':1,'lastFac':.01,'Nstep':8,'theta':1,'volt':3020,'linIter':200,'NewtTol':1e-6,'linTol':1e-8,'maxLoadIter':32,'np':'4~1~1'}]
   ls+=[{'testName':'pullin2D-lin-3050' ,'beamL':.1,'pois':.31,'Nthick':1,              'Nstep':8,'theta':1,'volt':3050,'linIter':200,'NewtTol':1e-6,'linTol':1e-8,'maxLoadIter':32,'np':'4~1~1'}]
   return ls


if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single

   #dataDict={'mainDir':mainDir,'pois':.31,'Nthick':2,'lastFac':.01,'Nstep':10,'theta':1,'volt':195,'linIter':10000,'NewtTol':1e-6,'linTol':1e-8}

   dataDict={'mainDir':mainDir,'beamL':.1,'pois':.31,'Nthick':1,'Nstep':10,'theta':1,'volt':2800,'linIter':200,'NewtTol':1e-6,'linTol':1e-8}

   single(dataDict)
   tipTrace()
   #animGif(dataDict['frameCount'])
   #animGif(3)











