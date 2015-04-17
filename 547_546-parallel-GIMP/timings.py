# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

# To run all tests and compare to standards, invoke this file as:
#    python timings.py
#
# To create a new standard for the tests invoke as:
#    python timings.py originate

import os,sys

#timingsDirList=['align','bar','cant3D-modes','fixed3D-tip','plasNorm','cant2D-modes','cant2D-body','fixed3D-surface','gen2-path','pullin2D'] # Add timing test directories to this list
timingsDirList=['align','fixed3D-surface']
#timingsDirList=['fixed3D-surface']

class timingMethods:
   nproc=1
   testTag='timing'
   @classmethod
   def getTestList(self,runDir):
      topDir=os.getcwd()
      try:
         os.chdir(runDir)             # cannot multi-thread this, or next
         sys.path.append(os.getcwd()) # allow importing of modules from here
         from run import timingTest
         testList=timingTest()    # returns candidate dictionaries
         sys.path.remove(os.getcwd()) # don't want to find this stuff again
         del sys.modules['run']       # must manually remove, else it stays in and other 'run' modules aren't recognized
         os.chdir(topDir)
         return testList
      except:
         print '    !! Could not getTestList() for',runDir
         traceback.print_exc(file=sys.stdout)
         os.chdir(topDir)
         return None
   @classmethod
   def compareDicts(self,src,dst,rf):
      ok=True
      keys=src.keys()
      keys.sort(key=str.lower)
      for key in keys:
         if key=='runDir':continue
         if key=='mainDir':continue
         if key=='wallTime':
            ct=float(dst[key])
            st=float(src[key])
            err=abs(ct-st)/st
            if err>0.01:
               print>>rf,"    Unexpected wall time; standard:",st,'  candidate:',ct
               ok=False
            continue
         if key=='estimated total wall time':continue
         if key=='estimated wall time remaining':continue
         if key in dst:
            if str(src[key])!=str(dst[key]):
               print>>rf,"    Failed comparison -- key:",key,"  standard:",src[key],'  candidate:',dst[key]
               ok=False
         else:
            print>>rf,"    Candidate missing key:",key,src[key]
            ok=False
      keys=dst.keys()
      keys.sort(key=str.lower)
      for key in keys:
         if key not in src:
            print>>rf,"    Standard missing key:",key,dst[key]
      if ok:print>>rf,'\t    ... Passed'
      else: print>>rf,'    !! Failed one or more comparisons !!'



mainDir=os.path.join(os.getcwd(),'0main') # Point to the main that is used to run the tests

if __name__ == "__main__":
   sys.path.append(mainDir) # allow importing of modules from mainDir
   import main
   timMeth=timingMethods
   if len(sys.argv)==2 and sys.argv[1]=="originate":
      main.makeStandards(timingsDirList,mainDir,timMeth)
   else:
      main.makeCandidates(timingsDirList,mainDir,timingMethods)


