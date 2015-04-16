// Philip Wallstedt 2004-2009
#ifndef SHAPE_H
#define SHAPE_H

#include "patch.h"
using namespace std;

// These integration and interpolation operators are common,
// regardless of what shape function is used.
struct shapeSC {
   virtual void integrate(const patch& pch, const vector<double>& pu, vector<double>& gu) = 0;
   virtual void integrate(const patch& pch, const vector<Vector2>& pu, vector<Vector2>& gu) = 0;
   virtual void divergence(const patch& pch, const vector<Matrix2>& pu, vector<Vector2>& gu) = 0;
   virtual void interpolate(const patch& pch, vector<Vector2>& pu, const vector<Vector2>& gu) = 0;
   virtual void gradient(const patch& pch, vector<Matrix2>& pu, const vector<Vector2>& gu) = 0;
   virtual void updateContribList(patch& pch) = 0;
   virtual int Nghost() = 0;
   virtual int Nsupport() = 0;
   virtual ~shapeSC() {}
};
shapeSC& makeShape(string s);

// The standard piecewise-linear shape functions used in MPM
struct tent {
   double Sx, Sy, Gx, Gy;
   int Nghost() {return 1;}
   int Nsupport() {return 4;}
   void uS(double& S, const double& x, const double& h) {S = 1. - abs(x) / h;}
   void uG(double& G, const double& x, const double& h) {G = -sgn(x) / h;}
   void updateSG(const patch& pch, const int p, const int m) {
      const double rx = pch.px[p].x - pch.gx[m].x;
      const double ry = pch.px[p].y - pch.gx[m].y;
      uS(Sx, rx, pch.dx);
      uS(Sy, ry, pch.dy);
      uG(Gx, rx, pch.dx);
      uG(Gy, ry, pch.dy);
   }
};

// The special adaptive spline shape functions for GIMP that are most accurate
struct GIMP {
   double Sx, Sy, Gx, Gy;
   int Nghost() {return 2;}
   int Nsupport() {return 9;}
   void uS(double& S, const double& x, const double& h, const double& l) {
      const double r = abs(x);
      if (r < l)  {S = 1. - (r * r + l * l) / (2.*h * l);         return;      }
      if (r < h - l) {S = 1. - r / h; return;}
      if (r < h + l) {S = (h + l - r) * (h + l - r) / (4.*h * l); return;}
      S = 0.;
   }
   void uG(double& G, const double& x, const double& h, const double& l) {
      const double r = abs(x);
      if (r < l) {G = -x / (h * l); return;}
      if (r < h - l) {G = -sgn(x) / h; return;}
      if (r < h + l) {G = (h + l - r) / (-2.*sgn(x) * h * l); return;}
      G = 0.;
   }
   void updateSG(const patch& pch, const int p, const int m) {
      const double rx = pch.px[p].x - pch.gx[m].x;
      const double ry = pch.px[p].y - pch.gx[m].y;
      const double hyx = pch.dy / pch.dx;
      // dimension factor for finding particle half-width
      // dfac=2 for 1D, dfac=4 for 2D, and dfac=8 for 3D
      const double dfac = 4.;
      const double lpxi = sqrt(pch.pV[p] / (dfac * pch.thick * hyx));
      const double lpyi = lpxi * hyx;
      uS(Sx, rx, pch.dx, lpxi * pch.pF[p].xx);
      uS(Sy, ry, pch.dy, lpyi * pch.pF[p].yy);
      uG(Gx, rx, pch.dx, lpxi * pch.pF[p].xx);
      uG(Gy, ry, pch.dy, lpyi * pch.pF[p].yy);
   }
};

// Each specialization of operators draws from these templated operators
template<typename S>struct operations: public shapeSC, public S {
   int Nghost() {return S::Nghost();}
   int Nsupport() {return S::Nsupport();}
   void integrate(const patch& pch, const vector<double>& pu, vector<double>& gu);
   void integrate(const patch& pch, const vector<Vector2>& pu, vector<Vector2>& gu);
   void divergence(const patch& pch, const vector<Matrix2>& pu, vector<Vector2>& gu);
   void interpolate(const patch& pch, vector<Vector2>& pu, const vector<Vector2>& gu);
   void gradient(const patch& pch, vector<Matrix2>& pu, const vector<Vector2>& gu);
   void updateContribList(patch& pch);
};

// Helper functions to create new operator instances
operations<tent>*newMPM();
operations<GIMP>*newGIMP();


#endif







