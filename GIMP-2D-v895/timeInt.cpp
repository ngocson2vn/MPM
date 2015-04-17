// Philip Wallstedt 2004-2009

#include "func.h"
using namespace std;

// There are several variations of the GIMP algorithm.
// We need to be able to try out new variations, but the best
// one so far is "cenDif".  However, "UVF" is needed for MPM.

void UVF::advance(patch& pch) {
   
   cerr << "\nUVF::advance " << pch.elaps << "\n";

   for (int i = 0; i < pch.Npart(); ++i) {
      pch.pJ[i] = getStress(pch, pch.pF[i], pch.pVS[i]);
      pch.pVS[i] *= pch.pV[i] * pch.pJ[i];
   }

   sh.integrate (pch, pch.pm, pch.gm);
   sh.integrate (pch, pch.pfe, pch.gfe);
   sh.divergence(pch, pch.pVS, pch.gfi);

   for (int i = 0; i < pch.Nnode(); ++i) pch.gm[i] += machTol;

   for (int i = 0; i < pch.Nnode(); ++i) pch.ga[i] = (pch.gfe[i] + pch.gfi[i]) / pch.gm[i];

   UVF::applyGridBC(pch);

   if (pch.incCount == 0) for (int i = 0; i < pch.Nnode(); ++i) pch.ga[i] *= .5;

   sh.interpolate(pch, pch.pvI, pch.ga);

   for (int i = 0; i < pch.Npart(); ++i) pch.pv[i] += pch.pvI[i] * pch.dt;

   for (int i = 0; i < pch.Npart(); ++i) pch.pw[i] = pch.pv[i] * pch.pm[i];

   sh.integrate (pch, pch.pw, pch.gw);

   for (int i = 0; i < pch.Nnode(); ++i) pch.gv[i] = pch.gw[i] / pch.gm[i];

   UVF::applyGridBC(pch);
   sh.interpolate(pch, pch.pxI, pch.gv);
   sh.gradient(pch, pch.pGv, pch.gv);

   for (int i = 0; i < pch.Npart(); ++i) pch.px[i] += pch.pxI[i] * pch.dt;

   for (int i = 0; i < pch.Npart(); ++i) pch.pF[i] += pch.pGv[i].inner(pch.pF[i]) * pch.dt;
}

void cenDif::advance(patch& pch) {
   // Find the stress and jacobian from the deformation gradient and user-supplied constitutive model (in the problem file)
   for (int i = 0; i < pch.Npart(); ++i) {
      pch.pJ[i] = getStress(pch, pch.pF[i], pch.pVS[i]);
      pch.pVS[i] *= pch.pV[i] * pch.pJ[i];
   }

   // Form momentum on each particle
   for (int i = 0; i < pch.Npart(); ++i) pch.pw[i] = pch.pv[i] * pch.pm[i];

   // Integrate mass to the grid, forming mass density
   sh.integrate (pch, pch.pm, pch.gm);

   // Integrate momentum to the grid, forming momentum density
   sh.integrate (pch, pch.pw, pch.gw);

   // Integrate external force to the grid
   sh.integrate (pch, pch.pfe, pch.gfe);

   // Integrate the divergence of stress to the grid, to get internal force
   sh.divergence(pch, pch.pVS, pch.gfi);

   // We'll divide by mass soon, so we don't want masses to be zero
   for (int i = 0; i < pch.Nnode(); ++i) pch.gm[i] += machTol;

   // Divide momentum by mass to get velocity
   for (int i = 0; i < pch.Nnode(); ++i) pch.gv[i] = pch.gw[i] / pch.gm[i];

   // Divide force by mass to get acceleration
   for (int i = 0; i < pch.Nnode(); ++i) pch.ga[i] = (pch.gfe[i] + pch.gfi[i]) / pch.gm[i];

   // Apply user-supplied boundary conditions to node velocities and accelerations
   cenDif::applyGridBC(pch);

   // Clever way of initializing the leap-frog central differencing
   if (pch.incCount == 0) for (int i = 0; i < pch.Nnode(); ++i) pch.ga[i] *= .5;

   // Advance the grid velocity
   for (int i = 0; i < pch.Nnode(); ++i) pch.gv[i] += pch.ga[i] * pch.dt;

   // Find the particle velocity increment
   sh.interpolate(pch, pch.pvI, pch.ga);

   // Find the particle position increment
   sh.interpolate(pch, pch.pxI, pch.gv);

   // Find the deformation gradient increment
   sh.gradient(pch, pch.pGv, pch.gv);

   // Advance the particle velocity
   for (int i = 0; i < pch.Npart(); ++i) pch.pv[i] += pch.pvI[i] * pch.dt;

   // Advance the particle position
   for (int i = 0; i < pch.Npart(); ++i) pch.px[i] += pch.pxI[i] * pch.dt;

   // Advance the deformation gradient
   for (int i = 0; i < pch.Npart(); ++i) pch.pF[i] += pch.pGv[i].inner(pch.pF[i]) * pch.dt;
}

timeIntSC& makeTimeInt(string s, shapeSC&sh, patch&) {
   if (s == "UVF") {
      return *new typename UVF::UVF(sh);
   }

   if (s == "CD" ) {
      return *new typename cenDif::cenDif(sh);
   }

   throw "makeTimeInt: no time integrator assigned";
}

