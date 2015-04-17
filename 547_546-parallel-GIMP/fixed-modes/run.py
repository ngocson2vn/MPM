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

def bicgstab():
   gf=open('gnucomm.xls','w')
   print>>gf,'set terminal png size 2000,1000 crop'
   print>>gf,"set key bottom left"
   print>>gf,"set grid"
   print>>gf,"set log y"
   print>>gf,"set xrange [0:20000]"
   print>>gf,"set xlabel 'Time'"
   print>>gf,"set ylabel 'Velocity'"
   print>>gf,"set output 'bicgstab.png'"
   print>>gf,"plot 'bicgstab.xls' using (abs($1)) title 'Kiter' with linespoints",
   print>>gf,",    'bicgstab.xls' using (abs($2)) title 'rho'   with linespoints",
   print>>gf,",    'bicgstab.xls' using (abs($3)) title 'beta'  with linespoints",
   print>>gf,",    'bicgstab.xls' using (abs($4)) title 'alpha' with linespoints",
   print>>gf,",    'bicgstab.xls' using (abs($5)) title 'dotTS' with linespoints",
   print>>gf,",    'bicgstab.xls' using (abs($6)) title 'dotTT' with linespoints",
   print>>gf,",    'bicgstab.xls' using (abs($7)) title 'omega' with linespoints",
   print>>gf,",    'bicgstab.xls' using (abs($8)) title 'ratio' with linespoints"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def regressionTest(): # should NOT specify mainDir
   ls=[]
   #ls+=[{'testName':'fixed-modes-lCG','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'linCG' ,'CFL':64,'load':1e-5,'mode':1,'linIter':400,'linIter':256,'linTol':1e-8}]
   #ls+=[{'testName':'fixed-modes-NCG','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'NewtCG','CFL':64,'load':1e-5,'mode':1,'linIter':400,'linIter':256,'linTol':1e-8}]
   #ls+=[{'testName':'fixed-modes-lGM','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'linGM' ,'CFL':64,'load':1e-5,'mode':1,'linIter':400,'linIter':256,'linTol':1e-8}]
   #ls+=[{'testName':'fixed-modes-NGM','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'dynGM' ,'CFL':64,'load':1e-5,'mode':1,'linIter':400,'linIter':256,'linTol':1e-8}]

   ls+=[{'testName':'fixed-modes-linGM','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'linGM' ,'CFL':64,'load':1e-5,'mode':1,'linIter':400,'linIter':256,'linTol':1e-8,'np':'2~1~1'}]
   ls+=[{'testName':'fixed-modes-dynCG','beamL':.1,'Nthick':2,'Nwide':1,'shape':'MPM','tInt':'NewtCG','CFL':128,'load':1e-5,'mode':1,'linIter':400,'NewtTol':1e-6,'linTol':1e-8,'np':'3~1~1'}]
   ls+=[{'testName':'fixed-modes-dynGM','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'dynGM','CFL':128,'load':1e-5,'mode':1,'linIter':256,'NewtTol':1e-6,'linTol':1e-8,'np':'4~1~1'}]
   ls+=[{'testName':'fixed-modes-linCG','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'linCG','CFL':64,'load':1e-5,'mode':1,'linIter':400,'linTol':1e-8,'np':'2~1~1'}]
   ls+=[{'testName':'fixed-modes-typical','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'momentum','CFL':.9,'load':1e-5,'mode':1,'np':'4~1~1'}]
   ls+=[{'testName':'fixed-modes-MPM-mom','beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'momentum','CFL':2.,'load':1e-5,'mode':1,'np':'3~1~1'}]
   ls+=[{'testName':'fixed-modes-MPM-CD' ,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'CD'      ,'CFL':.7,'load':1e-5,'mode':1,'np':'2~1~1'}]
   ls+=[{'testName':'fixed-modes-GIMP-mom','beamL':.1,'Nthick':1,'Nwide':1,'shape':'GIMP','tInt':'momentum','CFL':2.,'load':1e-5,'mode':1}]
   ls+=[{'testName':'fixed-modes-GIMP-CD' ,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'GIMP','tInt':'CD'      ,'CFL':.8,'load':1e-5,'mode':1}]
   ls+=[{'testName':'fixed-modes-spline-mom','beamL':.1,'Nthick':1,'Nwide':1,'shape':'spline','tInt':'momentum','CFL':2.,'load':1e-5,'mode':1}]
   ls+=[{'testName':'fixed-modes-spline-CD' ,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'spline','tInt':'CD'      ,'CFL':.8,'load':1e-5,'mode':1}]
   return ls


if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single

   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'dynBCG2','CFL':64,'load':1e-5,'mode':1,'linIter':400,'NewtTol':1e-6,'linTol':1e-8,'NewtLim':2}
   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':2,'Nwide':1,'shape':'MPM','tInt':'momGM','CFL':128,'load':1e-5,'mode':1,'linIter':400,'NewtTol':1e-6,'linTol':1e-8}
   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'momGM','CFL':128,'load':1e-5,'mode':1,'linIter':256,'NewtTol':1e-6,'linTol':1e-8}
   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'momentum','CFL':.9,'load':1e-5,'mode':1}
   dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'RKmom','CFL':2.,'load':1e-5,'mode':1}

   #dataDict={'mainDir':mainDir,'beamL':.1,'Nthick':1,'Nwide':1,'shape':'MPM','tInt':'NewtCG','CFL':64,'load':1e-5,'mode':1,'linIter':400,'linTol':1e-8}


   single(dataDict,False)
   #print 'wall time',dataDict['wallTime']
   tipTrace(dataDict)
   bicgstab()











