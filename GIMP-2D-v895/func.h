// Philip Wallstedt 2004-2009
#ifndef IO_H
#define IO_H
#include <sstream>
#include "patch.h"
#include "shape.h"
using namespace std;

// Convenience class for processing command line arguments
// If a match is found then the command line value is applied and
// a note is made to stdout.
// If no match is found then nothing happens, and there is no warning
class argMap: public map<string, string> {
   argMap::iterator it;
   string es;

public:
   argMap(int argc, char**argv);
   string operator[](const string& s) {
      it = find(s);
      if (it == end())return es;
      cerr << it->first << '=' << it->second << '\t';
      return it->second;
   }
};

// We let istringstream figure out the type based on where the argument will go
template<typename T> void operator>>(const string& s, T& t) {
   if (s.size() > 0) {
      istringstream i(s); i >> t;
   }
}

// Several algorithmic variations are used, but they all provide these functions
struct timeIntSC {
   virtual void advance(patch&) = 0;
   void applyGridBC(patch&); // problem specific
   double getStress(const patch&, const Matrix2&, Matrix2&); // problem specific
   virtual ~timeIntSC() {}
};

// io.cpp - utility functions for setting up arrangements of particles
void fillRectangleRegular(patch& pch, const Vector2& b, const Vector2& e, const double ppe, const int);
void fillAnnulusRegular(patch& pch, const Vector2& o, const double r, const double R, const double ppe, const int);
void fillAnnulusRadial(patch& pch, const Vector2& o, const double Ri, const double Ro, const int);

// timeInt.cpp - algorithmic variations
struct UVF : public timeIntSC {
   shapeSC& sh;
   
   void advance(patch& pch);

   UVF(shapeSC& s) : sh(s) {

   }
};

struct cenDif: public timeIntSC {shapeSC& sh; void advance(patch& pch); cenDif(shapeSC& s): sh(s) {}};
timeIntSC& makeTimeInt(string, shapeSC&, patch&);


#endif




