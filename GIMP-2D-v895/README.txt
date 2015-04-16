README for 2D GIMP code version 894
Philip Wallstedt - wallstedt@yahoo.com
17Feb2009

This software is released into the public domain; 
it may be used for any purpose, commercial or otherwise, free of charge.
No guarantee of accuracy or correctness is expressed or implied.
If you use this code, I would love to hear about it.
Write me at wallstedt@yahoo.com and tell me your experience
and I will try to answer any questions you might have.

System Recommendations:
- Python 2.5 or higher (not required, but without it you must provide your own command line arguments)
- Gnuplot 4 or higher (compiled with libgd) with the animated gif terminal (not required, but without it you must provide your own visualization)
- gcc 3.4.4 or higher (may work with Visual C++; the code only requires the C++ standard library)
- make (or compile patch.cpp, shape.cpp, io.cpp, timeInt.cpp, and a problem file yourself)

I have successfully run this code on cygwin (Windows 2000) and Linux and an older version ran on VC++.

Example of use where only gcc and make are available:
   > make impact
   > ./impact.exe Ncell=64 shape=GIMP tInt=CD CFL=.4 ppe=2 load=.1 p1=32

Example of use where Python, Gnuplot, and gcc are available:
   > python run.py  (this will compile, run, and plot the results in an animated .gif, which is best viewed in a web browser)

The meaning of the above command line parameters:
   Ncell   = characteristic number of cells, usually cells per dimension
   shape   = shape function used; may be piecewise linear "MPM", or "GIMP", or something else
   tInt    = time integration algorithm; recommend "CD"
   CFL     = Courant number; a measure of time step size; stable if less than .8 or so
   ppe     = Particles per dimension per cell; 2ppe gives 4 particles per cell, 3ppe = 9, etc.
   load    = Some characteristic loading parameter or force to the system
   p1      = particle visualization size; used with gnuplot

Other command line parameters may be used for different problems, and new parameters are easy to define with the argMap mechanism of func.h

Files in this package:
   Helpful utilities:
      makefile     : makes stuff
      run.py       : helps automate passing of command line parameters
   Core files:
      tensor.h     : defines vector and matrix operations
      patch.h      : the central repository for problem data
      patch.cpp    : patch constructor plus odds and ends
      shape.h      : shape functions and operator virtual functions
      shape.cpp    : instances of operator virtual functions
      func.h       : time integration virtual functions and utilities
      io.cpp       : functions to make various shapes of particles
      timeInt.cpp  : several time-stepping algorithms
   Problem files:
      impact.cpp   : shows off impact-modeling
      ring.cpp     : manufactured solution with moving surface
      align.cpp    : manufactured solution that should be 2nd order in GIMP

I made this code by reducing my research code down to essential elements while retaining place-holders
for things so that the framework and organization of the code is apparent.  The code
is designed for maximum flexibility so that new algorithmic variations are easy to
implement.  However, not all of the flexibility of this code is required for a useful code.
I've put comments in various files, mostly in patch.h and timeInt.cpp; the only
problem file with comments is align.cpp

The space() and time() functions in run.py allow measurements of convergence in codes
so correctness can be demonstrated, even if different computer architectures produce
results that are not numerically identical.  The align and ring problems have exact
solutions so they can be used for this purpose.  The impact problem shows off the
general ability of the GIMP algorithm.  The tear problem demonstrates a very simple
failure model.  The two(ball) problem appeared in the original paper for MPM.

The problem file is the highest level instruction for the code.
This version includes five problem files: impact, ring, align, tear, and two.
Each problem file defines the stress function,
   double timeIntSC::getStress(const patch&pch,const Matrix2&F,Matrix2&S)
the boundary condition function,
   void timeIntSC::applyGridBC(patch&pch)
a visualization output function like,
   void history(patch&pch,ostream&os,int&fc)
and may define an exact solution if one exists.  The problem file also
calls main and runs the general time stepping loop.

Each problem file assembles information from the core files
in roughly the following manner:
   1. #include "func.h" which includes "shape.h" and "patch.h"
   2. call main() and make a map of the command line arguments
   3. create the shape class shapeSC
   4. create one or more patches (just one for this code)
   5. create the time integration class timeIntSC
   6. Form a cloud of particles to represent a shape
   7. Iterate on the time step loop:
      1. Apply current body force to particles
      2. Advance forward by one time step
      3. Plot results
   8. Pass interesting results to stdout





