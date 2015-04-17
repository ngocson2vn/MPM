# Philip Wallstedt 2004-2009
from subprocess import *
import sys,time

def logLog(fn,Nser):
   gf=open('gnucomm.xls','w')
   xa='7' # "+xa+"
   ti='2' # "+ti+"
   print>>gf,'set terminal png size 640,480 crop font "Arial" 14'
   print>>gf,"set key top left"
   print>>gf,"set grid"
   print>>gf,"set log xy"
   print>>gf,"set xlabel 'N cells per dimension'"
   #print>>gf,"set xlabel 'CFL'"
   print>>gf,"set output 'L1.png'"
   print>>gf,"set ylabel 'L-1 error'"
   print>>gf,"plot '"+fn+"' index 0 using "+xa+":12 title "+ti+" with linespoints lw 4 ps 2",
   for n in range(1,Nser):
      print>>gf,", '"+fn+"' index "+str(n)+" using "+xa+":12 title "+ti+" with linespoints lw 4 ps 2",
   print>>gf
   print>>gf,"set output 'Li.png'"
   print>>gf,"set ylabel 'L-inf error'"
   print>>gf,"plot '"+fn+"' index 0 using "+xa+":14 title "+ti+" with linespoints lw 4 ps 2",
   for n in range(1,Nser):
      print>>gf,", '"+fn+"' index "+str(n)+" using "+xa+":14 title "+ti+" with linespoints lw 4 ps 2",
   gf.close()
   Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

# Measuring spatial convergence for the code by calling it many 
# times with gradually increasing mesh sizes.  Process handling
# is easy in Python and the whole thing is automated.  Gnuplot
# will put nice convergence plots in the current directory.
def space():
   case='align'
   CFL=['CFL=.4']
   if Popen(['make',case],shell=False).wait()!=0:sys.exit()
   serList=[]

   serList.append(['single=0','shape=GIMP'  ,'tInt=CD' ,'load=.1','ppe=2']+CFL)
#    serList.append(['single=0','shape=UGIMP'  ,'tInt=CD' ,'load=.1','ppe=2']+CFL)
#    serList.append(['single=0','shape=GIMP'  ,'tInt=USL' ,'load=.1','ppe=2']+CFL)
#    serList.append(['single=0','shape=GIMP'  ,'tInt=USF' ,'load=.1','ppe=2']+CFL)
#    serList.append(['single=0','shape=MPM'  ,'tInt=CD' ,'load=.1','ppe=2']+CFL)

   sf=open('space.keep','w')
   Nser=0
   for series in serList:
      Nser+=1
      print>>sf,'\n\n'
      try:
         header=[]
         for it in series:
            arg=it.split('=')
            if arg[0]=='tInt' or arg[0]=='shape':header+=[arg[1]]
            else:header+=[arg[0]]
         header+=['Ncell','Npart','Nnode','Nstep','wallTime','L1','L2','Linf']
         for r in header:print>>sf,r,'\t',
         print>>sf
         loN=10.
         hiN=100
         rN=1.7783
         N=loN
         while int(N)<=hiN:
            pars=series+['Ncell='+str(int(N))]
            prc=Popen(['./'+case+'.exe']+pars,shell=False,stdout=PIPE)
            prc.wait()
            vals=[]
            for it in pars:vals+=[it.split('=')[1]]
            try:vals+=prc.stdout.readlines()[0].split()
            except:print vals
            for r in vals:print>>sf,r,'\t',
            print>>sf
            sf.flush()
            N*=rN
      except KeyboardInterrupt:
         sf.flush()
         try:
            raw_input("Buffer flushed.  ENTER to continue with the next series or CTL-c again to stop gracefully:")
            continue
         except KeyboardInterrupt:break
   sf.close()
   logLog('space.keep',Nser)


# Measuring temporal convergence and plotting the result
def time():
   case='ring'
   Nc=['Ncell=56']
   if Popen(['make',case],shell=False).wait()!=0:sys.exit()
   serList=[]

   serList.append(['single=0','tInt=CD','shape=GIMP'  ,'load=.1','ppe=2']+Nc)
#    serList.append(['single=0','tInt=CD','shape=UGIMP'   ,'load=.1','ppe=2']+Nc)
#    serList.append(['single=0','tInt=USL','shape=GIMP'   ,'load=.1','ppe=2']+Nc)
#    serList.append(['single=0','tInt=USF','shape=GIMP'   ,'load=.1','ppe=2']+Nc)
#    serList.append(['single=0','tInt=CD','shape=MPM'   ,'load=.1','ppe=2']+Nc)

   sf=open('time.keep','w')
   Nser=0
   for series in serList:
      Nser+=1
      print>>sf,'\n\n'
      try:
         header=[]
         for it in series:
            arg=it.split('=')
            if arg[0]=='tInt' or arg[0]=='shape':header+=[arg[1]]
            else:header+=[arg[0]]
         header+=['CFL','Npart','Nnode','Nstep','wallTime','L1','L2','Linf']
         for r in header:print>>sf,r,'\t',
         print>>sf
         loN=.001
         hiN=1.
         rN=.05
         N=hiN
         while N>loN:
            pars=series+['CFL='+str(N)]
            prc=Popen(['./'+case+'.exe']+pars,shell=False,stdout=PIPE)
            prc.wait()
            vals=[]
            for it in pars:vals+=[it.split('=')[1]]
            try:vals+=prc.stdout.readlines()[0].split()
            except:print vals
            for r in vals:print>>sf,r,'\t',
            print>>sf
            sf.flush()
            N-=rN
      except KeyboardInterrupt:
         sf.flush()
         try:
            raw_input("Buffer flushed.  ENTER to continue with the next series or CTL-c again to stop gracefully:")
            continue
         except KeyboardInterrupt:break
   sf.close()
   logLog('time.keep',Nser)

def ring():
   args=['Ncell=16','shape=GIMP','tInt=CD' ,'CFL=.4','ppe=2','load=.2','p1=32']
   # run make; stop if an error is found
   if Popen(['make','ring'],shell=False).wait()!=0:sys.exit()
   # run the program and pick up stdout
   prc=Popen(['./ring.exe']+args,shell=False,stdout=PIPE)
   prc.wait()
   # read the number of frames from stdout
   ni=int(prc.stdout.readlines()[0].rstrip())
   # make an animated gif of the results
   print "Requested Frames:",ni,
   gf=open('gnucomm.xls','w')
   print>>gf,"set term gif animate opt delay 10 size 800,600 crop background '#7D96A4' font 'tahoma' 14"
   print>>gf,"set output 'ring.gif'"
   print>>gf,"set view map"
   print>>gf,"set cbrange [-8000:6000]"
   print>>gf,"set xrange [0:1.2]"
   print>>gf,"set yrange [0:1.2]"
   print>>gf,"set size square"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,'set title "'+str(i)+' / '+str(ni)+'"'
      print>>gf,"splot 'history.xls' index "+str(i)+" using 1:2:3:($4) with points pt 5 ps var lt palette"
   gf.close()
   Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def impact():
   args=['Ncell=32','shape=GIMP','tInt=CD' ,'CFL=.4','ppe=2','load=15','p1=24']
   # run make; stop if an error is found
   if Popen(['make','impact'],shell=False).wait()!=0:sys.exit()
   # run the program and pick up stdout
   prc=Popen(['./impact.exe']+args,shell=False,stdout=PIPE)
   prc.wait()
   # read the number of frames from stdout
   ni=int(prc.stdout.readlines()[0].rstrip())
   # make an animated gif of the results
   print "Requested Frames:",ni,
   gf=open('gnucomm.xls','w')
   print>>gf,"set term gif animate opt delay 10 size 800,600 crop background '#7D96A4' font 'tahoma' 14"
   print>>gf,"set output 'impact.gif'"
   print>>gf,"set view map"
   print>>gf,"set cbrange [.5:1.5]"
   print>>gf,"set xrange [0:2]"
   print>>gf,"set yrange [0:1]"
   print>>gf,"set size ratio .5"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,'set title "'+str(i)+' / '+str(ni)+'"'
      print>>gf,"splot 'history.xls' index "+str(i)+" using 1:2:3:($4) with points pt 5 ps var lt palette"
   gf.close()
   Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def align():
   args=['Ncell=16','shape=GIMP','tInt=CD','CFL=.4','ppe=2','load=.1','p1=24']
   # run make; stop if an error is found
   if Popen(['make','align'],shell=False).wait()!=0:sys.exit()
   # run the program and pick up stdout
   prc=Popen(['./align.exe']+args,shell=False,stdout=PIPE)
   prc.wait()
   # read the number of frames from stdout
   ni=int(prc.stdout.readlines()[0].rstrip())
   # make an animated gif of the results
   print "Requested Frames:",ni,
   gf=open('gnucomm.xls','w')
   print>>gf,"set term gif animate opt delay 10 size 800,600 crop background '#7D96A4' font 'tahoma' 14"
   print>>gf,"set output 'align.gif'"
   print>>gf,"set view map"
   print>>gf,"set cbrange [.5:1.5]"
   print>>gf,"set xrange [0:1]"
   print>>gf,"set yrange [0:1]"
   print>>gf,"set size square"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,'set title "'+str(i)+' / '+str(ni)+'"'
      print>>gf,"splot 'history.xls' index "+str(i)+" using 1:2:3:($4) with points pt 5 ps var lt palette"
   gf.close()
   Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def tear():
   args=['Ncell=32','shape=GIMP','tInt=CD' ,'CFL=.4','ppe=2','load=15','p1=24','ult=1.5e6']
   # run make; stop if an error is found
   if Popen(['make','tear'],shell=False).wait()!=0:sys.exit()
   # run the program and pick up stdout
   prc=Popen(['./tear.exe']+args,shell=False,stdout=PIPE)
   prc.wait()
   # read the number of frames from stdout
   ni=int(prc.stdout.readlines()[0].rstrip())
   # make an animated gif of the results
   print "Requested Frames:",ni,
   gf=open('gnucomm.xls','w')
   print>>gf,"set term gif animate opt delay 10 size 800,600 crop background '#7D96A4' font 'tahoma' 14"
   print>>gf,"set output 'tear.gif'"
   print>>gf,"set view map"
   #print>>gf,"set cbrange [.5:1.5]"
   print>>gf,"set xrange [0:2]"
   print>>gf,"set yrange [0:1]"
   print>>gf,"set size ratio .5"
   print>>gf,'set palette rgbformulae 33,13,10'
   print>>gf,"set nokey"
   print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,'set title "'+str(i)+' / '+str(ni)+'"'
      print>>gf,"splot 'history.xls' index "+str(i)+" using 1:2:3:($4) with points pt 5 ps var lt palette"
   gf.close()
   Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

def two():
   args=['Ncell=32','shape=MPM','tInt=UVF','CFL=.4','ppe=2','load=.1','p1=8'] # MPM with UVF
   #args=['Ncell=32','shape=GIMP','tInt=CD','CFL=.4','ppe=2','load=.1','p1=8'] # GIMP with CD
   # run make; stop if an error is found
   if Popen(['make','two'],shell=False).wait()!=0:sys.exit()
   # run the program and pick up stdout
   prc=Popen(['./two.exe']+args,shell=False,stdout=PIPE)
   prc.wait()
   # read the number of frames from stdout
   ni=int(prc.stdout.readlines()[0].rstrip())
   # make an animated gif of the results
   print "Requested Frames:",ni,
   gf=open('gnucomm.xls','w')
   #print>>gf,"set term gif animate opt delay 60 size 800,600 crop background '#7D96A4' font 'arial' 14"
   # print>>gf,"set term pngcairo size 800,600 crop background '#7D96A4' font 'arial' 14"
   # print>>gf,"set output 'two.gif'"
   # print>>gf,"set view map"
   # print>>gf,"set cbrange [0:50]"
   # print>>gf,"set xrange [0:1]"
   # print>>gf,"set yrange [0:1]"
   # print>>gf,"set size square"
   # print>>gf,'set palette rgbformulae 33,13,10'
   # print>>gf,"set nokey"
   # print>>gf,"set border 0"
   for i in xrange(ni):
      print>>gf,"set term png size 800,600 crop background '#7D96A4'"
      print>>gf,"set output 'output/two%s.png'" % i
      print>>gf,"set view map"
      print>>gf,"set cbrange [0:50]"
      print>>gf,"set xrange [0:1]"
      print>>gf,"set yrange [0:1]"
      print>>gf,"set size square"
      print>>gf,'set palette rgbformulae 33,13,10'
      print>>gf,"set nokey"
      print>>gf,"set border 0"
      print>>gf,'set title "' + str(i)+ ' / ' + str(ni) + '"'
      print>>gf,"splot 'history.xls' index " + str(i) + " using 1:2:3:($4) with points pt 5 ps var palette"
   gf.close()
   Popen(['gnuplot','gnucomm.xls'],shell=False).wait()

#space()
#time()

#ring()
#impact()
#align()
#tear()
two()





