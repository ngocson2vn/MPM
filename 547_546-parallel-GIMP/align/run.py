# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

import os,sys,subprocess

def logLog(fn,Nser):
   gf=open('gnucomm.xls','w')
   xa='3' # "+xa+"
   ti='1' # "+ti+"
   print>>gf,'set terminal png size 640,480 crop font "Arial" 14'
   print>>gf,"set key top left"
   print>>gf,"set grid"
   print>>gf,"set log xy"
   print>>gf,"set xlabel 'N cells per dimension'"
   #print>>gf,"set xlabel 'CFL'"
   print>>gf,"set output 'L1.png'"
   print>>gf,"set ylabel 'L-1 error'"
   print>>gf,"plot '"+fn+"' index 0 using "+xa+":8 title "+ti+" with linespoints lw 4 ps 2",
   for n in range(1,Nser):
      print>>gf,", '"+fn+"' index "+str(n)+" using "+xa+":8 title "+ti+" with linespoints lw 4 ps 2",
   print>>gf
   print>>gf,"set output 'Li.png'"
   print>>gf,"set ylabel 'L-inf error'"
   print>>gf,"plot '"+fn+"' index 0 using "+xa+":10 title "+ti+" with linespoints lw 4 ps 2",
   for n in range(1,Nser):
      print>>gf,", '"+fn+"' index "+str(n)+" using "+xa+":10 title "+ti+" with linespoints lw 4 ps 2",
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

# Measuring spatial convergence for the code by calling it many times with gradually increasing mesh sizes.
def space(mainDir):
   sf=open('space.keep','w')
   print>>sf,'\n\n'
   header=['tInt','shape','Ncell','Npart','Nnode','incCount','wallTime','L1norm','L2norm','Linfnorm']
   for r in header:print>>sf,r,'\t',
   print>>sf
   loN=10.
   hiN=90.
   rN=1.7783
   N=loN
   while int(N)<=hiN:
      try:
         #outDict={'mainDir':mainDir,'Ncell':int(N),'single':0,'shape':'spline','tInt':'CD','load':.1,'CFL':.7}
         outDict={'mainDir':mainDir,'Ncell':int(N),'single':0,'shape':'spline','tInt':'CD','load':.1,'CFL':.7,'np':'2~2~1'}
         single(outDict,False)
         for r in header:print>>sf,outDict[r],'\t',
         print>>sf
         sf.flush()
         N*=rN
      except KeyboardInterrupt:
         sf.flush()
         break
   sf.close()
   logLog('space.keep',1)

def animGif(outDict):
   ni=int(outDict['frameCount'])
   print "Requested Frames:",ni,' ',
   gf=open('gnucomm.xls','w')
   print>>gf,"set term gif animate opt delay 10 size 900,900 crop x000000"
   print>>gf,"set output 'movie.gif'"
   print>>gf,"set view map"
   print>>gf,"set border linecolor rgb 'cyan'"
   #print>>gf,"set cbrange [0:100]"
   print>>gf,"set xrange [0:1]"
   print>>gf,"set yrange [0:1]"
   print>>gf,"set size square"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,'set title "'+str(i)+' / '+str(ni)+'" tc rgb "cyan"'
      print>>gf,"splot 'history.xls' index "+str(i)+" using 1:2:4 with points pt 5 ps 1 lt palette"
   gf.close()
   subprocess.Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def timingTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'align-GIMP-CD','Ncell':23,'shape':'GIMP','tInt':'CD' ,'CFL':.9,'load':.1}]
   ls+=[{'testName':'align-MPM-mom','Ncell':27,'shape':'MPM','tInt':'momentum','CFL':.9,'load':.1}]
   return ls

def regressionTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'align-GIMP-CD','Ncell':8,'shape':'GIMP','tInt':'CD'      ,'CFL':.9,'load':.1}]
   ls+=[{'testName':'align-MPM-mom','Ncell':8,'shape':'MPM' ,'tInt':'momentum','CFL':.9,'load':.1,'np':'2~2~2'}]
   return ls



if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single
   dataDict={'mainDir':mainDir,'Ncell':8,'shape':'GIMP','tInt':'CD','CFL':.5,'load':.1,'np':'2~2~2'}
   #dataDict={'mainDir':mainDir,'Ncell':5,'shape':'MPM','tInt':'momentum' ,'CFL':.9,'load':.1}
   #single(dataDict,False)
   #animGif(dataDict)
   space(mainDir)
   #logLog('space.keep',1)



