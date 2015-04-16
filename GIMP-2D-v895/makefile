# Philip Wallstedt 2004-2009

CXX = g++
CXXFLAGS = -O3 -ffast-math -Wall -Wextra -Wno-uninitialized -I/usr/include/c++/4.8.3/i686-redhat-linux
COMMONFILES = patch.o shape.o io.o timeInt.o

impact: $(COMMONFILES) impact.o
	$(CXX) $(CXXFLAGS) -o impact.exe $(COMMONFILES) impact.o
	strip impact.exe

ring: $(COMMONFILES) ring.o
	$(CXX) $(CXXFLAGS) -o ring.exe $(COMMONFILES) ring.o
	strip ring.exe

align: $(COMMONFILES) align.o
	$(CXX) $(CXXFLAGS) -o align.exe $(COMMONFILES) align.o
	strip align.exe

tear: $(COMMONFILES) tear.o
	$(CXX) $(CXXFLAGS) -o tear.exe $(COMMONFILES) tear.o
	strip tear.exe

two: $(COMMONFILES) two.o
	$(CXX) $(CXXFLAGS) -o two.exe $(COMMONFILES) two.o
	strip two.exe

clean:
	rm -f *.o *.exe *.stackdump

version:
	rm -f *.exe *.stackdump *.xls *.png *.pdf *.gif *~

patch.o:		tensor.h patch.h
shape.o:		tensor.h patch.h shape.h
io.o:       tensor.h patch.h func.h
timeInt.o:  tensor.h patch.h func.h
impact.o:   tensor.h patch.h func.h
ring.o:     tensor.h patch.h func.h
align.o:    tensor.h patch.h func.h
tear.o:     tensor.h patch.h func.h
two.o:      tensor.h patch.h func.h






