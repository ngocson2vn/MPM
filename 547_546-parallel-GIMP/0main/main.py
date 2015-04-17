# Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
# This software is covered by the license and liability disclaimer in "license.txt"

import sys,time,math,os,cPickle,subprocess,string,traceback,os.path,Tkinter,threading
import multiprocessing as mp
from Queue import Empty as emptyQ

NcompileProc=8

def runOnly(di):
   args=['mpirun']
   if('np' in di):
      nx,ny,nz=str(di['np']).split('~')
      args+=['-np',str(int(nx)*int(ny)*int(nz))]
   else:args+=['-np','1']
   args+=[os.path.join(di['runDir'],'run.exe')]
   for d in di.items():
      if d[0]=='runDir':continue
      if d[0]=='mainDir':continue
      args+=[d[0]+'='+str(d[1])]
   argStr=string.join(args)
   prc=subprocess.Popen(argStr,shell=True,cwd=di['runDir'])
   prc.wait()
   if 'testName' in di:repName=di['testName']+'.report'
   else:repName='run.report'
   reportPath=os.path.join(di['runDir'],repName)
   outList=open(reportPath,'r').readlines()
   for line in outList:
      tpl=line.strip('\n').partition('=')
      if tpl[1]=='':di[tpl[0]]=None
      else:di[tpl[0]]=tpl[2]
   try:
      if 'testName' in di:fitName=di['testName']+'.fit'
      else:fitName='run.fit'
      reportPath=os.path.join(di['runDir'],fitName)
      outList=open(reportPath,'r').readlines()
      for line in outList:
         tpl=line.strip('\n').partition('=')
         if tpl[1]=='':di[tpl[0]]=None
         else:di[tpl[0]]=tpl[2]
   except:pass

def runQueue(jobs,numberOfThreads):
   def worker():
      while True:
         try:
            probDict=jobs.get_nowait()
            runOnly(probDict)        # here is where the work gets done
            completed.put(probDict)  # collect completed jobs
            jobs.task_done()
         except emptyQ:return        # this is normal behavior
         except:                     # this is abnormal
            jobs.task_done()         # we MUST end the failed task or the Queue will hang
            print 'Job skipped',probDict
            traceback.print_exc(file=sys.stdout)
   completed=mp.Queue()
   lenjobs=jobs.qsize()
   for i in range(numberOfThreads):
      t=mp.Process(target=worker)
      t.start()
   jobs.join()                       # wait until all queue items have called task_done
   series=[]
   while completed.qsize() > 0:series.append(completed.get_nowait())
   return series

def oneList(listOfDicts,mainDir,np):
   runDir=os.getcwd()
   if subprocess.Popen('scons / -j '+str(NcompileProc)+' runDir='+runDir,shell=True,cwd=mainDir).wait()!=0:sys.exit()
   idx=0
   jobs=Queue.Queue()
   for dict in listOfDicts:
      dict['testName']=str(idx)
      print dict
      dict['runDir']=runDir
      jobs.put(dict)
      idx+=1
   doneList=runQueue(jobs,np) # do all the work - this will take some time
   return sorted(doneList,key=lambda dict:dict['testName'])

def runGUI():
   class tagVal:
      def __init__(self,row):
         self.tag = Tkinter.StringVar()
         self.tagLab = Tkinter.Label(textvariable=self.tag,anchor="e",fg="lightblue",bg="black",font=("Helvetica",12),height=1,width=32)
         self.tagLab.grid(column=0,row=row)
         self.val = Tkinter.StringVar()
         self.valLab = Tkinter.Label(textvariable=self.val,anchor="w",fg="lightblue",bg="black",font=("Helvetica",12),height=1,width=32)
         self.valLab.grid(column=1,row=row)
      def set(self,tag,val):
         self.tag.set(tag)
         self.val.set(val)
   go=True
   tvList=[]
   if os.path.exists('run.report'):os.remove('run.report')
   prog=Tkinter.Tk()
   prog.title("Working Directory: "+str(os.getcwd()))
   prog.grid()
   tvList+=[tagVal(1)]
   tvList[0].set("Starting ","Progress Indicator")
   while go:
      lineList=[]
      for chances in range(4):
         if os.path.exists('run.report'):
            progFile=open('run.report','r')
            lineList=progFile.readlines()
            progFile.close()
            if len(lineList)>0:break
         time.sleep(1)
      if len(tvList)!=len(lineList):
         tvList=[]
         for i in range(len(lineList)):tvList+=[tagVal(len(tvList)+1)]
      for i in range(len(lineList)):
         tpl=lineList[i].strip('\n').partition('=')
         tvList[i].set(tpl[0]+': ',tpl[2])
         #if tpl[0]=='wallTime':go=False
         if tpl[0]=='wallTime':return
      prog.update()
      time.sleep(1)
   #prog.mainloop()

def single(di,showGUI=True):
   di['runDir']=os.getcwd()
   if subprocess.Popen('scons / -j '+str(NcompileProc)+' runDir='+di['runDir'],shell=True,cwd=di['mainDir']).wait()!=0:sys.exit()
   if showGUI:gui=threading.Thread(target=runGUI)
   if showGUI:gui.start()
   runOnly(di)
   cPickle.dump(di,open("dataDict",'w'))
   #if 'wallTime' not in di.keys():gui.join()
   if showGUI:return gui

def runCandidates(dirList,mainDir,testKind):
   topDir=os.getcwd()
   candDicts=[]
   for dir in dirList: # visit every directory in the list
      runDir=os.path.join(topDir,dir)
      candList=testKind.getTestList(runDir)
      for dict in candList:
         dict['np']='1~1~1' # quick hack to run all regressions in serial
         dict['runDir']=runDir
         dict['mainDir']=mainDir
      candDicts+=candList # add them to the general list
      if subprocess.Popen('scons / -j '+str(NcompileProc)+' runDir='+runDir,shell=True,cwd=mainDir).wait()!=0:sys.exit()
   jobs=mp.JoinableQueue()
   for dict in candDicts:
      print dict['testName'],dict['runDir']
      jobs.put(dict)
   doneList=runQueue(jobs,testKind.nproc) # do all the work - this will take some time
   doneDict={}
   for dict in doneList:
      doneDict[dict['testName']]=dict # a dictionary of dictionaries, indexed by testName
   return doneDict

def makeStandards(dirList,mainDir,testKind):
   topDir=os.getcwd()
   stanDict=runCandidates(dirList,mainDir,testKind) # run all tests
   import traceback
   for dir in dirList: # visit every directory in the list
      try:
         runDir=os.path.join(topDir,dir)
         localList=testKind.getTestList(runDir)
         localSave=[]
         print "   Making "+testKind.testTag+" standards in",dir
         for test in localList:
            print '     ',test['testName']
            localSave.append(stanDict[test['testName']]) # collect the tests invoked from this dir
         savePath=os.path.join(runDir,testKind.testTag+'Tests')
         cPickle.dump(localSave,open(savePath,'w')) # write standard file in dir
      except:
         print '    !! Failed, exception thrown, original standard could not be created'
         traceback.print_exc(file=sys.stdout)
         continue

def makeCandidates(dirList,mainDir,testKind):
# First, lets verify that we have standards for the tests we plan to run.
   topDir=os.getcwd()
   stanDict={}
   import traceback
   OK=True
   for dir in dirList:                      # visit every directory in the list
      try:
         runDir=os.path.join(topDir,dir)
         savePath=os.path.join(runDir,testKind.testTag+'Tests')
         localStanList=cPickle.load(file(savePath,'r'))
         for dict in localStanList:
            stanDict[dict['testName']]=dict # add to dictionary of dictionaries, indexed by testName
         localCandList=testKind.getTestList(runDir)
         localCandDict={}
         for dict in localCandList:localCandDict[dict['testName']]=dict # so we can locate by name
         for dict in localStanList:
            if dict['testName'] in localCandDict:pass
            else:
               print 'Candidate',dict['testName'],'is missing from run.py: '+testKind.testTag+'Test()'
               OK=False
         for dict in localCandList:
            if dict['testName'] in stanDict:pass
            else:
               print 'Standard',dict['testName'],'is missing from',savePath
               OK=False
      except:
         traceback.print_exc(file=sys.stdout)
         print '    !! Failed, exception thrown in makeCandidates (early loop) from',dir
   if not OK:
      raw_input("Hit ENTER to proceed with "+testKind.testTag+" tests.")
# Then, run the tests
   candDict=runCandidates(dirList,mainDir,testKind) # run all tests - most work done here
# Finally, compare the tests with the standards
   rf=open(testKind.testTag+'TestReport','w')
   for stan in stanDict.values():          # compare all tests
      print>>rf,testKind.testTag+' test for',stan['testName'],
      try:
         testKind.compareDicts(stan,candDict[stan['testName']],rf)
      except:
         traceback.print_exc(file=sys.stdout)
         print>>rf,'    !! Failed, exception thrown, test could not be performed'
   rf.close()
   rf=open(testKind.testTag+'TestReport','r')
   print rf.read()

