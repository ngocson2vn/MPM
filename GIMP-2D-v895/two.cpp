
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
      double CFL = .2;                   am["CFL"  ] >> CFL;
      am["single"] >> single;
      shapeSC& shp = makeShape(shpstr);
      patch pch(Nc, Nc, 0., 0., 1., 1., shp.Nghost(), 1.);
      timeIntSC&ti = makeTimeInt(tmistr, shp, pch);
      pch.ppe = 2.;                      am["ppe"  ] >> pch.ppe;
      pch.load = .1;                     am["load" ] >> pch.load;
      pch.Ymod = 1000.;
      pch.dens = 1000.;
      pch.vwav = sqrt(pch.Ymod / pch.dens);
      pch.pois = .3;
      pch.damp = .0;
      pch.partSize = 32.;                am["p1"   ] >> pch.partSize;
      pch.dt = min(pch.dx, pch.dy) * CFL / sqrt(pch.Ymod / pch.dens);
      fillAnnulusRegular(pch, Vector2(.25, .25), 0., .2, pch.ppe, shp.Nsupport());
      fillAnnulusRegular(pch, Vector2(.75, .75), 0., .2, pch.ppe, shp.Nsupport());
      pch.elaps = .0 / pch.vwav;

      for (int i = 0; i < pch.Npart(); i += 1) {
         pch.px[i] = pch.pX[i];
         pch.pv[i] = (pch.px[i].y > 1. - pch.px[i].x ? -1. : 1.) * Vector2(pch.load, pch.load);
         pch.pF[i] = I2();
      }

      shp.updateContribList(pch);

      // main loop
      ofstream hf("history.xls");
      int frameCount = 0;
      double Li = 0., L2 = 0., L1 = 0.;
      
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
      catch (const char*s) {cerr << "Defined Exception:" << s << endl; L1 = L2 = Li = 1.;}
      catch (const std::exception&error) {cerr << "Standard Exception:" << error.what() << endl;}
      catch (...) {cerr << "Unknown Exception" << endl;}
      cerr << '\n';

      // post process
      cerr << "last step: " << pch.incCount << endl;

      if (single) {
         cerr << "Wall time:" << double(clock() - t) / CPS << tab << "L-inf:" << Li << endl;
         cout << frameCount;
      } else {
         cerr << "Wall time:" << double(clock() - t) / CPS << tab << "L-inf:" << Li << endl;
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