// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include <fstream>
#include "main.h"
#include "constit.h"
using namespace std;

struct uniaxialXupdatedElastic:public constitutiveSC{
   vector<mat3>stress0;
   vector<realT>volume0;
   const realT Gmod,Kmod,Dens,Pois,twoThirds;
   const partArray<mat3>&velGrad;
   partArray<mat3>&defGrad;
   partArray<realT>&volume;
   partArray<mat3>&stress;
   void update(const realT&);
   void revert();
   void save();
   realT waveSpeed()const{return sqrt((Kmod+Gmod*7./3.)/Dens);}
   uniaxialXupdatedElastic(const realT Y,
                 const realT P,
                 const realT D,
                 const partArray<mat3>&vg,
                 partArray<mat3>&dg,
                 partArray<realT>&vl,
                 partArray<mat3>&st):
                     Gmod(Y/(2.*(1.+P))),
                     Kmod(Y/(3.*(1.-2.*P))),
                     Dens(D),
                     Pois(P),
                     twoThirds(2./3.),
                     velGrad(vg),
                     defGrad(dg),
                     volume(vl),
                     stress(st){}
};

// could optimize revert so arrays are swapped rather than copied
void uniaxialXupdatedElastic::revert(){
   for(int i=0;i<Npart();++i)stress[i]=stress0[i];
   for(int i=0;i<Npart();++i)volume[i]=volume0[i];
}

// could optimize save so resize is only called once
void uniaxialXupdatedElastic::save(){
   stress0.resize(Npart());
   volume0.resize(Npart());
   for(int i=0;i<Npart();++i)stress0[i]=stress[i];
   for(int i=0;i<Npart();++i)volume0[i]=volume[i];
}

void uniaxialXupdatedElastic::update(const realT&delt){
   for(int i=0;i<Npart();++i){
      defGrad[i]+=velGrad[i].inner(defGrad[i])*delt;
      mat3 strainInc(.5*(velGrad[i]+trans(velGrad[i]))*delt);
      strainInc.yy=-Pois*strainInc.xx;
      strainInc.zz=-Pois*strainInc.xx;
      const realT evol=trace(strainInc);
      volume[i]+=evol*volume[i];
      stress[i]+=2.*Gmod*strainInc+I3((Kmod-twoThirds*Gmod)*evol);
   }
}

struct myPatch:public patch{
   const shapeSC&shp;
   const realT beamL,beamH,beamB,Ymod,dens,pois,load;
   updatedElastic*cstPtr;
   ofstream hf;
   int frameCount,probe1Q,probe2Q,probe3Q,probe4Q;
   bool afterStep();
   void gridVelocityBC(nodeArray<vec3>&,const realT&)const;
   myPatch(constitPtr&,const shapeSC&,const realT,const realT,const realT);
};

void myPatch::gridVelocityBC(nodeArray<vec3>&v,const realT&)const{
   for(int i=0;i<Nnode();++i){
      if(shp.curPosN[i].x<machEps*beamL)v[i].x=0.;
   }
}

bool myPatch::afterStep(){
   //hf<<elapsedTime<<tab
     //<<curPosP[probe1Q].x-refPosP[probe1Q].x<<tab
     //<<curPosP[probe2Q].x-refPosP[probe2Q].x<<tab
     //<<curPosP[probe3Q].x-refPosP[probe3Q].x<<tab
     //<<curPosP[probe4Q].x-refPosP[probe4Q].x<<nwl;
   hf<<elapsedTime<<tab
     <<stressP[probe1Q].xx<<tab
     <<stressP[probe2Q].xx<<tab
     <<stressP[probe3Q].xx<<tab
     <<stressP[probe4Q].xx<<nwl;
   return elapsedTime<finalTime;
   //return false;
}

myPatch::myPatch(constitPtr&cst,
                 const shapeSC&s,
                 const realT bl,
                 const realT bh,
                 const realT bb):
                 patch(s),
                 shp(s),
                 beamL(bl),
                 beamH(bh),
                 beamB(bb),
                 Ymod(comLineArg("Ymod",1e3)),
                 dens(comLineArg("dens",10.)),
                 pois(comLineArg("pois",.25)),
                 load(comLineArg("load",1.)){
   frameCount=0;
   cst=new uniaxialXupdatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   //cst=new updatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   elapsedTime=0.;
   finalTime=4.*beamL/sqrt(Ymod/dens);
   report.param("finalTime",finalTime);
   vec3 ppd=vec3(2.,1.,1.);
   fillBox(*this,shp,vec3(0.,0.,0.),vec3(beamL,beamH,beamB),ppd);
   for(int i=0;i<Npart();++i){
      velP[i]=load;
      massP[i]=refVolP[i]*dens;
   }
   probe1Q=partProbe(*this,.5*vec3( .5*beamL,beamH,beamB));
   probe2Q=partProbe(*this,.5*vec3(    beamL,beamH,beamB));
   probe3Q=partProbe(*this,.5*vec3(1.5*beamL,beamH,beamB));
   probe4Q=partProbe(*this,.5*vec3(2. *beamL,beamH,beamB));
   hf.open("history.xls");  hf.precision(14);
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   const int Nlong=comLineArg("Nlong",10.);
   const vec3 lo(0.,0.,0.);
   const vec3 hi(2.,1.,1.);
   makeShape(shp,2*Nlong,1,1,lo,hi,comLineArg("shape",string("MPM")));
   pch=new myPatch(cst,*shp,1.,1.,1.);
   makeTimeInt(tmi,pch,cst,shp,comLineArg("tInt",string("momentum")));
}












