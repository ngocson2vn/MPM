// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include "constit.h"
using namespace std;

// could optimize revert so arrays are swapped rather than copied
void neoHookean::revert(){
   for(indexT i=0;i<Npart();++i)defGrad[i]=defGrad0[i];
   for(indexT i=0;i<Npart();++i)stress[i]=stress0[i];
   for(indexT i=0;i<Npart();++i)volume[i]=volume0[i];
}

// could optimize save so resize is only called once
void neoHookean::save(){
   defGrad0.resize(Npart());
   stress0.resize(Npart());
   volume0.resize(Npart());
   for(indexT i=0;i<Npart();++i)defGrad0[i]=defGrad[i];
   for(indexT i=0;i<Npart();++i)stress0[i]=stress[i];
   for(indexT i=0;i<Npart();++i)volume0[i]=volume[i];
}

void neoHookean::update(const realT&delt){
   for(indexT i=0;i<Npart();++i){
      defGrad[i]+=velGrad[i].inner(defGrad[i])*delt;
      const realT J=det(defGrad[i]);
      //volume[i]+=trace(velGrad[i])*volume[i]*delt; // first order with align3D
      volume[i]=refVol[i]*J;                         // second order with align3D
      stress[i]=I3(lam*log(J)/J)+(mu/J)*(defGrad[i].inner(trans(defGrad[i]))-I3());
   }
}


// could optimize revert so arrays are swapped rather than copied
void updatedElastic::revert(){
   for(indexT i=0;i<Npart();++i)defGrad[i]=defGrad0[i];
   for(indexT i=0;i<Npart();++i)stress[i]=stress0[i];
   for(indexT i=0;i<Npart();++i)volume[i]=volume0[i];
}

// could optimize save so resize is only called once
void updatedElastic::save(){
   defGrad0.resize(Npart());
   stress0.resize(Npart());
   volume0.resize(Npart());
   for(indexT i=0;i<Npart();++i)defGrad0[i]=defGrad[i];
   for(indexT i=0;i<Npart();++i)stress0[i]=stress[i];
   for(indexT i=0;i<Npart();++i)volume0[i]=volume[i];
}

void updatedElastic::update(const realT&delt){
   for(indexT i=0;i<Npart();++i){
      defGrad[i]+=velGrad[i].inner(defGrad[i])*delt; // Nanson works for all
      const mat3 strainInc(.5*(velGrad[i]+trans(velGrad[i]))*delt);
      //defGrad[i]+=strainInc; // Nanson works with shear problem, but not with fixed-fixed beam
      const realT evol=trace(strainInc);
      volume[i]+=evol*volume[i];
      stress[i]+=2.*shearModG*strainInc+I3((bulkModK-twoThirds*shearModG)*evol);
   }
}


// could optimize revert so arrays are swapped rather than copied
void J2plasticLinearIsoKin::revert(){
   for(indexT i=0;i<Npart();++i)defGrad[i]=defGrad0[i];
   for(indexT i=0;i<Npart();++i)stress[i]=stress0[i];
   for(indexT i=0;i<Npart();++i)volume[i]=volume0[i];
   for(indexT i=0;i<Npart();++i)strain[i]=strain0[i];
   for(indexT i=0;i<Npart();++i)plasticStrain[i]=plasticStrain0[i];
   for(indexT i=0;i<Npart();++i)backStress[i]=backStress0[i];
   for(indexT i=0;i<Npart();++i)internalAlpha[i]=internalAlpha0[i];
}

// could optimize save so resize is only called once
void J2plasticLinearIsoKin::save(){
   defGrad0.resize(Npart());
   stress0.resize(Npart());
   volume0.resize(Npart());
   for(indexT i=0;i<Npart();++i)defGrad0[i]=defGrad[i];
   for(indexT i=0;i<Npart();++i)stress0[i]=stress[i];
   for(indexT i=0;i<Npart();++i)volume0[i]=volume[i];
   for(indexT i=0;i<Npart();++i)strain0[i]=strain[i];
   for(indexT i=0;i<Npart();++i)plasticStrain0[i]=plasticStrain[i];
   for(indexT i=0;i<Npart();++i)backStress0[i]=backStress[i];
   for(indexT i=0;i<Npart();++i)internalAlpha0[i]=internalAlpha[i];
}

void J2plasticLinearIsoKin::update(const realT&delt){
   for(indexT i=0;i<Npart();++i){
      defGrad[i]+=velGrad[i].inner(defGrad[i])*delt;
      const mat3 strainInc(.5*(velGrad[i]+trans(velGrad[i]))*delt);
      const realT evol=trace(strainInc);
      volume[i]+=volume[i]*evol;
      strain[i]+=strainInc;
      const mat3 devStrainInc(dev(strainInc));
      const mat3 dilStressInc(I3(bulkModK*evol));
      const mat3 trialStress(dev(stress[i])+2.*shearModG*devStrainInc);
      const realT hardLaw=yieldStress+isoHardMod*internalAlpha[i];
      const realT ffYieldCond=Frob(trialStress)-sqrtTwoThirds*hardLaw;
      if(ffYieldCond<=0.){
         stress[i]+=2.*shearModG*devStrainInc+dilStressInc;
      }
      else{
         const mat3 unitStress(trialStress/Frob(trialStress));
         const realT gammaInc=ffYieldCond/(2.*(shearModG+oneThird*(isoHardMod+kinHardMod)));
         const mat3 plasticStrainInc(gammaInc*unitStress);
         internalAlpha[i]+=sqrtTwoThirds*gammaInc;
         backStress[i]+=(sqrtTwoThirds*kinHardMod*gammaInc)*unitStress;
         plasticStrain[i]+=plasticStrainInc;
         stress[i]+=2.*shearModG*(devStrainInc-plasticStrainInc)+dilStressInc;
      }
   }
}




