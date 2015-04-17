
=============== README for MPM-GIMP ================

This is a parallelized implementation of the Material Point Method (MPM) [Sulsky et al 1994, 1995] and Generalized Interpolation Material Point method (GIMP) [Bardenhagen, Kober 2004].  The code is written in C++, parallelized with MPI, built with Scons and run from Python.  The parallelization work is recent and has not been tested in large scale computations.  Output is to text files; no binary format has yet been implemented.

There are three options for shape functions: piecewise linear, second degree splines, and cpGIMP; however, only the piecewise linear functions are currently parallelized; the others will follow.

Two explicit time integration algorithms are included: the USL and momentum formulations.  Several parallel matrix-free implicit solvers are available including linear and nonlinear versions of Conjugate Gradient and GMRES, in dynamic and quasi-static forms.

The code is bundled with nearly thirty regression tests to ensure that it remains high quality even as additional features are added.  The code design has been as careful and standardized as possible, making use of the modularization and abstraction mechanisms of C++ and keeping public and private data separate.  The core code is bundled together in a static library, then linked with a separate problem file that provides initial and boundary conditions and forcing functions.  This design avoids the need to "hardwire" special conditions into the core code, and gives maximum flexibility to describe the problem at hand.

System Requirements:
- C++98 (or better) standards-compliant compiler (g++, msvc++, intel, etc.)
- MPI 1.1 or better
- Python 2.4 or better
- Scons 1.2 or better

System Recommendations:
- Gnuplot 4 or higher (compiled with libgd) with the animated gif terminal; several of the provided problems generate output plots.
- Paraview for visualization (using the legacy .vtk text file format)

If Scons is outdated or unavailable, it can easily be installed in your user space as follows:
1. Download the scons tarball scons.xx.tar.gz.
2. Unpack it in the current directory with tar -xvf scons.xx.tar.gz.
3. cd into the scons.xx directory.
4. Run the following command which will create /bin and /lib directories in your home directory (if they are not already present) and install scons there:
> python setup.py install --prefix=$HOME
5. Add ~/bin (created above) to your path.

With all software installed the entire code may be compiled and tested against existing regression results by going to the base directory of the code and typing
> python regress.py

Individual problems are run (and compiled) by going into the problem directory and typing
> python run.py

The BSD license is used because I want the code to be compatible with GPL'd projects and also with proprietary projects.  Comments and advice regarding the license for this code are welcome.

If you use this code, I would love to hear about it.  Write me at wallstedt@yahoo.com and tell me your experience and I will help you get the code working and offer advice on how to use it for your application.











