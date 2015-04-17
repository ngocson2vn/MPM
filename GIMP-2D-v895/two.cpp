
// Philip Wallstedt 2004-2009
// This problem is taken from Sulsky, Chen and Schreyer
// "A particle method for history-dependent materials" CMAME 1994.

#include <fstream>
#include <ctime>
#include <cmath>
#include "func.h"
using namespace std;

static bool single = true;
static double globalLiError = 0.;
static double globalL2Error = 0.;
static double globalL1Error = 0.;

double planeStrainNeoHookean(const patch& pch, const Matrix2& F, Matrix2& S) {
   const double v = pch.pois;
   const double l = pch.Ymod * v / ((1. + v) * (1. - 2.*v));
   const double m = .5 * pch.Ymod / (1. + v);
   const double Ja = det(F);
   
   if (Ja < machTol) {
      cerr << F << tab << Ja << endl; 
      throw "negative Jacobian";
   }

   S = I2(l * log(Ja) / Ja);
   S += m / Ja * (F.inner(trans(F)) - I2());

   return Ja;
}

double timeIntSC::getStress(const patch& pch, const Matrix2& F, Matrix2& S) {
   return planeStrainNeoHookean(pch, F, S);
}

double neoHookEnergy(const patch& pch, const Matrix2& F) {
   const double v = pch.pois;
   const double lam = pch.Ymod * v / ((1. + v) * (1. - 2.*v));
   const double mu = .5 * pch.Ymod / (1. + v);
   const double J = det(F);
   const double lnJ = log(J);
   const Matrix2 C = trans(F).inner(F);

   return .5 * lam * lnJ * lnJ - mu * lnJ + .5 * mu * (trace(C) - 2.);
}

void history(patch& pch, ostream& os, int& fc) { // animGif, radialGif, etc
   if (single && (pch.incCount % 8 == 0)) {
      fc += 1;
      os << '#' << fc << '\n';

      for (int i = 0; i < pch.Npart(); i += 1) {
         //os<<pch.px[i]<<'\t'<<pch.pJ[i]<<'\t'<<pch.partSize/double(pch.I)<<'\n';
         //cerr << "px[i]: " << pch.px[i] << "\n";
         os << pch.px[i] << '\t' << neoHookEnergy(pch, pch.pF[i]) << '\t' << pch.partSize / double(pch.I) << '\n';
      }

      os << "\n\n";
   }
}

void timeIntSC::applyGridBC(patch&) {}

int main(int argc, char**argv) {
   try {
      const double CPS = double(CLOCKS_PER_SEC);
      clock_t t = clock();

      // assemble command line parameter map
      argMap am(argc, argv);

      int Nc = 20;                       am["Ncell"] >> Nc;
      string shpstr = "GIMP";            am["shape"] >> shpstr;
      string tmistr = "cen";             am["tInt" ] >> tmistr;
      double CFL = 0.2;                  am["CFL"  ] >> CFL;
      am["single"] >> single;

      // create the shape class shapeSC
      shapeSC& shp = makeShape(shpstr);

      // create one patch
      patch pch(Nc, Nc, 0.0, 0.0, 1.0, 1.0, shp.Nghost(), 1.0);
      
      // create the time integration class timeIntSC
      timeIntSC& ti = makeTimeInt(tmistr, shp, pch);

      pch.ppe  = 2.0;                    am["ppe"  ] >> pch.ppe;
      pch.load = 0.1;                    am["load" ] >> pch.load;
      pch.Ymod = 1000.0;
      pch.dens = 1000.0;
      pch.vwav = sqrt(pch.Ymod / pch.dens);
      pch.pois = 0.3;
      pch.damp = 0.0;
      pch.partSize = 32.0;               am["p1"   ] >> pch.partSize;
      pch.dt = min(pch.dx, pch.dy) * CFL / sqrt(pch.Ymod / pch.dens);

      // Form a cloud of particles to represent two disks

      // disk 1
      fillAnnulusRegular(pch, Vector2(0.25, 0.25), 0.0, 0.2, pch.ppe, shp.Nsupport());

      // disk 2
      fillAnnulusRegular(pch, Vector2(0.75, 0.75), 0.0, 0.2, pch.ppe, shp.Nsupport());

      pch.elaps = 0.0 / pch.vwav;

      // Set initial velocity and initial deformation gradient
      for (int i = 0; i < pch.Npart(); i += 1) {
         // Initial position
         pch.px[i] = pch.pX[i];

         // Initial volocity
         pch.pv[i] = (pch.px[i].y > 1.0 - pch.px[i].x ? -1.0 : 1.0) * Vector2(pch.load, pch.load);

         // Initial deformation gradient
         pch.pF[i] = I2();
      }

      shp.updateContribList(pch);

      // main loop
      ofstream hf("history.xls");
      int frameCount = 0;
      double Li = 0.0, L2 = 0.0, L1 = 0.0;
      
      cerr << '\n';

      try {
         while (true) {
            ti.advance(pch);
            pch.elaps += pch.dt;
            pch.incCount += 1;
            history(pch, hf, frameCount);
            L1 = (globalL1Error > L1 ? globalL1Error : L1);
            L2 = (globalL2Error > L2 ? globalL2Error : L2);
            Li = (globalLiError > Li ? globalLiError : Li);
            
            for (int i = 0; i < pch.Npart(); ++i) {
               if (!pch.inRegion(pch.px[i])) {
                  throw "special out";
               }
            }

            shp.updateContribList(pch);
            cerr << '~';
         }
      }
      catch (const char*s) {cerr << "\nDefined Exception: " << s << endl; L1 = L2 = Li = 1.;}
      catch (const std::exception&error) {cerr << "\nStandard Exception: " << error.what() << endl;}
      catch (...) {cerr << "\nUnknown Exception" << endl;}
      cerr << '\n';

      // post process
      cerr << "last step: " << pch.incCount << endl;

      if (single) {
         cerr << "Wall time: " << double(clock() - t) / CPS << tab << "L-inf: " << Li << endl;
         cout << frameCount;
      } else {
         cerr << "Wall time: " << double(clock() - t) / CPS << tab << "L-inf: " << Li << endl;
         cout << pch.Npart() << '\t'
              << pch.I*pch.J << '\t'
              << pch.incCount << '\t'
              << double(clock() - t) / CPS << '\t'
              << L1 << '\t'
              << L2 << '\t'
              << Li;
      }

   } catch (...) {cerr << "Non-loop Exception" << endl;}
}