// Philip Wallstedt 2004-2009

#include "shape.h"
using namespace std;

template<typename S>
void setWeightGrad(S& s, patch& pch, const int p, const int cnt, const int k) {
   s.updateSG(pch, p, k);
   const double w = s.Sx * s.Sy;
   const double x = s.Gx * s.Sy;
   const double y = s.Gy * s.Sx;
   pch.pCon[p].portionArray[cnt].setPortion(k, w, x, y);
}

// Here we load the contrib list for each particle - a list of
// which nodes this particle contributes to, and the weights for each
template<typename S>
void operations<S>::updateContribList(patch& pch) {
   S s;
   assert(S::Nsupport() == 4 || S::Nsupport() == 9);
   const int& I = pch.I;
   if (S::Nsupport() == 4) {
      for (int p = 0; p < pch.Npart(); p++) {
         const int k = pch.inCell(pch.px[p]);
         setWeightGrad(s, pch, p, 0, k    );
         setWeightGrad(s, pch, p, 1, k + 1  );
         setWeightGrad(s, pch, p, 2, k + I  );
         setWeightGrad(s, pch, p, 3, k + I + 1);
         pch.pCon[p].Npor = 4;
      }
   }
   if (S::Nsupport() == 9) {
      for (int p = 0; p < pch.Npart(); p++) {
         const int k = pch.inCell9(pch.px[p]);
         setWeightGrad(s, pch, p, 0, k      );
         setWeightGrad(s, pch, p, 1, k + 1    );
         setWeightGrad(s, pch, p, 2, k + 2    );
         setWeightGrad(s, pch, p, 3, k + I    );
         setWeightGrad(s, pch, p, 4, k + I + 1  );
         setWeightGrad(s, pch, p, 5, k + I + 2  );
         setWeightGrad(s, pch, p, 6, k + I + I  );
         setWeightGrad(s, pch, p, 7, k + I + I + 1);
         setWeightGrad(s, pch, p, 8, k + I + I + 2);
         pch.pCon[p].Npor = 9;
      }
   }
}

// Below are defined the main operators that form the core of
// the GIMP algorithm - the integration and interpolation tasks.
template<typename S>
void operations<S>::integrate(const patch& pch, const vector<double>& pu, vector<double>& gu) {
   for (unsigned g = 0; g < gu.size(); g++)gu[g] = 0.;
   for (unsigned p = 0; p < pu.size(); p++) {
      const partContribs& pcon = pch.pCon[p];
      for (int k = 0; k < pcon.Npor; k += 1) {
         const partContribs::portion& por = pcon[k];
         gu[por.idx] += pu[p] * por.weight;
      }
   }
}

template<typename S>
void operations<S>::integrate(const patch& pch, const vector<Vector2>& pu, vector<Vector2>& gu) {
   for (unsigned g = 0; g < gu.size(); g++)gu[g] = 0.;
   for (unsigned p = 0; p < pu.size(); p++) {
      const partContribs& pcon = pch.pCon[p];
      const Vector2& puR = pu[p];
      for (int k = 0; k < pcon.Npor; k += 1) {
         const partContribs::portion& por = pcon[k];
         Vector2& guR = gu[por.idx];
         guR += puR * por.weight;
      }
   }
}

template<typename S>
void operations<S>::interpolate(const patch& pch, vector<Vector2>& pu, const vector<Vector2>& gu) {
   for (unsigned p = 0; p < pu.size(); p++)pu[p] = 0.;
   for (unsigned p = 0; p < pu.size(); p++) {
      const partContribs& pcon = pch.pCon[p];
      Vector2& puR = pu[p];
      for (int k = 0; k < pcon.Npor; k += 1) {
         const partContribs::portion& por = pcon[k];
         const Vector2& guR = gu[por.idx];
         puR.x += guR.x * por.weight;
         puR.y += guR.y * por.weight;
      }
   }
}

template<typename S>
void operations<S>::gradient(const patch& pch, vector<Matrix2>& pu, const vector<Vector2>& gu) {
   for (unsigned p = 0; p < pu.size(); p++)pu[p] = 0.;
   for (unsigned p = 0; p < pu.size(); p++) {
      const partContribs& pcon = pch.pCon[p];
      Matrix2& puR = pu[p];
      for (int k = 0; k < pcon.Npor; k += 1) {
         const partContribs::portion& por = pcon[k];
         const Vector2& guR = gu[por.idx];
         puR.xx += guR.x * por.gradx;
         puR.xy += guR.x * por.grady;
         puR.yx += guR.y * por.gradx;
         puR.yy += guR.y * por.grady;
      }
   }
}

template<typename S>
void operations<S>::divergence(const patch& pch, const vector<Matrix2>& pu, vector<Vector2>& gu) {
   for (unsigned g = 0; g < gu.size(); g++)gu[g] = 0.;
   for (unsigned p = 0; p < pu.size(); p++) {
      const partContribs& pcon = pch.pCon[p];
      const Matrix2& puR = pu[p];
      for (int k = 0; k < pcon.Npor; k += 1) {
         const partContribs::portion& por = pcon[k];
         Vector2& guR = gu[por.idx];
         guR.x -= puR.xx * por.gradx + puR.yx * por.grady;
         guR.y -= puR.xy * por.gradx + puR.yy * por.grady;
      }
   }
}


// Definitions of helper functions
operations<tent>* newMPM() {return new operations<tent>;}
operations<GIMP>* newGIMP() {return new operations<GIMP>;}


// At the beginning of the program we select which shape function
// we will use and send out its pointer to whoever wants it.
// Hopefully the virtual function penalty is only paid once.
shapeSC& makeShape(string s) { // want single v-table lookup
   shapeSC*t = NULL;
   if (s == "MPM") {t = newMPM();}
   if (s == "GIMP") {t = newGIMP();}
   if (t == NULL)throw"bad shape string";
   return *t;
}




