# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

import os,sys,subprocess

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
   Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def regressionTest(): # should NOT specify mainDir
   ls=[]
   ls+=[{'testName':'pressBox-CG','Ncell':5,'Nstep':8,'pres':6e9,'tInt':'quasiCG','linIter':600,'NewtTol':1e-6,'linTol':1e-8}]
   ls+=[{'testName':'pressBox-GM','Ncell':5,'Nstep':4,'pres':6e9,'tInt':'quasi'  ,'linIter':300,'NewtTol':1e-6,'linTol':1e-8}]
   ls+=[{'testName':'pressBox-BCG-6','Ncell':5,'Nstep':16,'pres':6e9,'tInt':'quasiBCG','linIter':8,'NewtTol':1e-100,'linTol':1e-100}]
   ls+=[{'testName':'pressBox-BCG-12','Ncell':5,'Nstep':16,'pres':12e9,'tInt':'quasiBCG','linIter':8,'NewtTol':1e-100,'linTol':1e-100}]
   return ls



if __name__ == "__main__":
   topDir=os.path.split(os.getcwd())[0]
   mainDir=os.path.join(topDir,'0main') # Define the path to the main program files
   sys.path.append(mainDir)             # allow importing of modules from mainDir
   from main import single
   #dataDict={'mainDir':mainDir,'Ncell':8,'shape':'GIMP','tInt':'CD' ,'CFL':.9,'load':.1}
   dataDict={'mainDir':mainDir,'Ncell':5,'Nstep':16,'pres':6e9,'tInt':'quasiBCG','linIter':32,'NewtTol':1e-6,'linTol':1e-100}
   #dataDict={'mainDir':mainDir,'Ncell':8,'Nstep':20,'pres':10e9,'tInt':'quasi','linIter':400,'NewtTol':1e-8,'linTol':1e-6}
   single(dataDict,False)
   #for el in outDict:print el,'=',outDict[el]
   #animGif(outDict)
   #space(mainDir)
   #logLog('space.keep',1)



