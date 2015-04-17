// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef CONSTIT_H
#define CONSTIT_H
#include "main.h"
using namespace std;

struct neoHookean:public constitutiveSC{
   vector<mat3>defGrad0,stress0;
   vector<realT>volume0;
   const realT Ymod,Pois,Dens,lam,mu;
   const partArray<mat3>&velGrad;
   partArray<mat3>&defGrad;
   partArray<realT>&volume;
   partArray<realT>&refVol;
   partArray<mat3>&stress;
   void update(const realT&);
   void revert();
   void save();
   realT waveSpeed()const{return sqrt((lam+3.*mu)/Dens);}
   neoHookean(const realT Y,
              const realT P,
              const realT D,
              const partArray<mat3>&vg,
              partArray<mat3>&dg,
              partArray<realT>&vl,
              partArray<realT>&rv,
              partArray<mat3>&st):
                     Ymod(Y),
                     Pois(P),
                     Dens(D),
                     lam(Ymod*Pois/((1.+Pois)*(1.-2.*Pois))),
                     mu(.5*Ymod/(1.+Pois)),
                     velGrad(vg),
                     defGrad(dg),
                     volume(vl),
                     refVol(rv),
                     stress(st){
      //if(Npart()<1)throw stop("neoHookean:no particles to initialize!");
   }
};

struct updatedElastic:public constitutiveSC{
   vector<mat3>stress0,defGrad0;
   vector<realT>volume0;
   const realT twoThirds,bulkModK,shearModG,density;
   const partArray<mat3>&velGrad;
   partArray<mat3>&defGrad;
   partArray<realT>&volume;
   partArray<mat3>&stress;
   void update(const realT&);
   void revert();
   void save();
   realT waveSpeed()const{return sqrt((bulkModK+shearModG*7./3.)/density);}
   updatedElastic(const realT Y,
                  const realT P,
                  const realT D,
                  const partArray<mat3>&vg,
                  partArray<mat3>&dg,
                  partArray<realT>&vl,
                  partArray<mat3>&st):
                     twoThirds(2./3.),
                     bulkModK(Y/(3.*(1.-2.*P))),
                     shearModG(Y/(2.*(1.+P))),
                     density(D),
                     velGrad(vg),
                     defGrad(dg),
                     volume(vl),
                     stress(st){
   }
};

struct J2plasticLinearIsoKin:public constitutiveSC{
   partArray<mat3>strain,plasticStrain,backStress;
   partArray<mat3>strain0,plasticStrain0,backStress0;
   partArray<realT>internalAlpha;
   partArray<realT>internalAlpha0;
   vector<mat3>stress0,defGrad0;
   vector<realT>volume0;
   const realT oneThird,sqrtTwoThirds,density,bulkModK,shearModG,yieldStress,isoHardMod,kinHardMod;
   const partArray<mat3>&velGrad;
   partArray<mat3>&defGrad;
   partArray<realT>&volume;
   partArray<mat3>&stress;
   void update(const realT&);
   void revert();
   void save();
   realT waveSpeed()const{return sqrt((bulkModK+shearModG*7./3.)/density);}
   J2plasticLinearIsoKin(const realT Y,
                         const realT P,
                         const realT D,
                         const realT c,
                         const realT d,
                         const realT e,
                         const partArray<mat3>&vg,
                         partArray<mat3>&dg,
                         partArray<realT>&vl,
                         partArray<mat3>&st):
                            oneThird(1./3.),
                            sqrtTwoThirds(sqrt(2./3.)),
                            density(D),
                            bulkModK(Y/(3.*(1.-2.*P))),
                            shearModG(Y/(2.*(1.+P))),
                            yieldStress(c),
                            isoHardMod(d),
                            kinHardMod(e),
                            velGrad(vg),
                            defGrad(dg),
                            volume(vl),
                            stress(st){
      for(indexT i=0;i<Npart();++i){
         internalAlpha[0]=0.;
      }
   }
};

/////////////  Schreyer  ////////////////
struct anisotropicElastic:public constitutiveSC{
   orth6 EE;
   mat3 stress,strainInc;
   anisotropicElastic(const vec3&Ymod,const vec3&shearG,const mat3&pois){
      mat3 F;
      for(indexT i=0;i<3;++i){
         for(indexT j=0;j<3;++j)F(i,j)=(i==j?1.:-1.)*pois(i,j)/Ymod[j];
      }
      EE=orth6(inv(F),2.*shearG);
   }
   void update(const realT&){
      stress+=EE.inner(strainInc);
   }
   realT waveSpeed()const{return 0.;}
};

struct isoElastic:public constitutiveSC{
   const realT twoThirds,bulkModK,shearModG;
   mat3 stress,strainInc;
   isoElastic(const realT&Y,const realT&P):twoThirds(2./3.),
                                           bulkModK(Y/(3.*(1.-2.*P))),
                                           shearModG(Y/(2.*(1.+P))){}
   mat3 EgivenS(const mat3&S){
      assert(S.xy==S.yx);
      assert(S.yz==S.zy);
      assert(S.xz==S.zx);
      const realT A=1./(2.*shearModG);
      const realT B=(1./(9.*bulkModK)-1./(6.*shearModG))*trace(S);
      return S*A+I3()*B;
   }
   void update(const realT&){
      const realT evol=trace(strainInc);
      stress+=2.*shearModG*strainInc+I3((bulkModK-twoThirds*shearModG)*evol);
   }
   realT waveSpeed()const{return 0.;}
};


#endif

