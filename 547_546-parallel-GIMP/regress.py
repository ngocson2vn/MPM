# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

# To run all tests and compare to standards, invoke this file as:
#    python regress.py
#
# To create a new standard for the tests invoke as:
#    python regress.py originate

import os,sys

regressDirList=['align','fixed-modes','fixed-uniform','gen2-path','plasNorm','pressBox','pullin2D','resStress','three'] # Add regression test directories to this list
#regressDirList=['align']

# replace oldString with newString recursively to *
# rpl -R -x'.cpp' -x'.h' -x'.py' 'oldString' 'newString' *

skipKeys=['np',
          'runDir',
          'mainDir',
          'wallTime',
          'estimated total wall time',
          'estimated wall time remaining',
          'connectivity table size',
          'inner connectivity table size',
          'outer connectivity table size']

class regressMethods:
   nproc=1
   nproc=8 # quick hack for serial regressions
   testTag='regression'
   @classmethod
   def getTestList(self,runDir):
      topDir=os.getcwd()
      try:
         os.chdir(runDir)             # cannot multi-thread this, or next
         sys.path.append(os.getcwd()) # allow importing of modules from here
         from run import regressionTest
         testList=regressionTest()    # returns candidate dictionaries
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
         if key in skipKeys:continue
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
         if key in skipKeys:continue
         if key not in src:
            print>>rf,"    Standard missing key:",key,dst[key]
      if ok:print>>rf,'\t    ... Passed'
      else: print>>rf,'    !! Failed one or more comparisons !!'



mainDir=os.path.join(os.getcwd(),'0main') # Point to the main that is used to run the tests

if __name__ == "__main__":
   sys.path.append(mainDir) # allow importing of modules from mainDir
   import main
   regMeth=regressMethods
   if len(sys.argv)==2 and sys.argv[1]=="originate":
      main.makeStandards(regressDirList,mainDir,regMeth)
   else:
      main.makeCandidates(regressDirList,mainDir,regressMethods)


