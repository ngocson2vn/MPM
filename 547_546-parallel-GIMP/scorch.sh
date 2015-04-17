#!/bin/bash
find * -type f -name "*.exe" -exec rm -fv {} \;
find * -type f -name "*.xls" -exec rm -fv {} \;
find * -type f -name "*.png" -exec rm -fv {} \;
find * -type f -name "*.pdf" -exec rm -fv {} \;
find * -type f -name "*.gif" -exec rm -fv {} \;
find * -type f -name "*.pyc" -exec rm -fv {} \;
find * -type f -name "*.out" -exec rm -fv {} \;
find * -type f -name "*.fit" -exec rm -fv {} \;
find * -type f -name "*.vtk" -exec rm -fv {} \;
find * -type f -name "*.o"   -exec rm -fv {} \;
find * -type f -name "*.a"   -exec rm -fv {} \;
find * -type f -name "*~"    -exec rm -fv {} \;
find * -type f -name "*.report" -exec rm -fv {} \;
find * -type f -name "*.dblite" -exec rm -fv {} \;
