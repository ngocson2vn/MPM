// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

struct myPatch:public patch{
   mutable partArray<vec3>surfaceP,posIncP;
   timeIntPtr&tmi;
   const shapeSC&shp;
   const realT beamL,beamH,beamB;
   vec3 ppd;
   realT Ymod,dens,pois;
   mutable realT currVolt;
   realT finalVolt,gap0;
   ofstream hf,tf,vf,nf;
   int frameCount;
   realT defGivenPres(const realT fperA){                    // fixed-fixed beam with uniform load
      const realT I=beamB*beamH*beamH*beamH/12.;
      const realT w=-fperA*beamB;                            // load per length
      const realT x=.5*beamL;                                // midpoint of beam
      return w*x*x*(beamL-x)*(beamL-x)/(24.*Ymod*I);
   }
   realT presGivenFracDef(const realT fracDef){              // fixed-fixed beam with uniform load
      const realT I=beamB*beamH*beamH*beamH/12.;
      const realT deflect=fracDef*beamH/beamB;               // find load required to deflect by this amount
      const realT x=.5*beamL;                                // midpoint of beam
      return 24.*Ymod*I*deflect/(x*x*(beamL-x)*(beamL-x));
   }
   realT presGivenGap(const realT gap_mm,const realT volt)const{  // fixed-fixed beam with electro-static load
      const realT gap=1e-3*gap_mm;                           // (meter)
      const realT permittivity=8.854187817e-12;              // (Farad/meter)
      const realT fperA=permittivity*volt*volt/(2.*gap*gap); // (Newton/meter^2)
      return 1e-6*fperA;                                     // (Newton/mm)
   }
   bool afterStep();
   void history();
   void gridVelocityBC(nodeArray<vec3>&,const realT&)const;
   void makeExternalForces(nodeArray<vec3>&,const realT&)const;
   myPatch(timeIntPtr&,constitPtr&,shapePtr&,const realT,const realT,const realT);
   ~myPatch();
};

void myPatch::gridVelocityBC(nodeArray<vec3>&v,const realT&)const{
   for(indexT i=0;i<Nnode();++i){
      if(shp.gridPosN(i).x<machEps*beamL || shp.gridPosN(i).x>oneEps*beamL)v[i]=0.;
   }
}

void myPatch::makeExternalForces(nodeArray<vec3>&gravityN,const realT&delt)const{
   const realT loadFactor=(tmi->elapsedTime() - tmi->timeStep() + delt) / tmi->finalTime;
   currVolt=finalVolt*loadFactor;
   report.param("loadFactor",loadFactor);
   shp.interpolate(posIncP,velN);
   realT mingap=huge;
   const vec3 snapNorm(0,-1,0);
   for(indexT i=0;i<Npart();++i)if(surfNormP[i].y<machEps){
      const realT refArea=surfAreaP[i].inner(snapNorm);
      const vec3 trialPosP=curPosP[i]+posIncP[i] * tmi->timeStep();
      const realT gap=gap0-abs(trialPosP.y-refPosP[i].y);
      if(gap<mingap)mingap=gap;
      if(gap<gap0*machEps)throw stop("Pull-in detected; exiting early");
      const realT pressure=presGivenGap(gap,currVolt);
      surfaceP[i]=(pressure*det(defGradP[i])*refArea*snapNorm.inner(inv(defGradP[i])));
   }
   shp.integrate(surfaceP,gravityN);
   for(indexT i=0;i<Nnode();++i)gravityN[i]/=massN[i];
   //report.param("current gap",mingap);
}

void myPatch::history(){
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1){
      tf<<currVolt<<tab<<gap0-abs(curPosP[i].y-refPosP[i].y)<<nwl;
   }

   frameCount+=1;
   hf<<'#'<<frameCount<<nwl;
   for(indexT i=0;i<Npart();++i)if(abs(refPosP[i].z)<.5*shp.cell.z)hf<<curPosP[i]<<tab<<vonMises(stressP[i])<<nwl;
   hf<<nwl<<endl;

   for(indexT i=0;i<Nnode();++i)nf<<shp.gridPosN(i)<<tab<<0<<nwl;
   nf<<nwl<<endl;

   //for(indexT i=0;i<surfList.size();++i)if(surfList[i].norm.y==-1.){
      //const partSurface&ps=surfList[i];
      //if(abs(refPosP[ps.p].z)>.5*shp.cell.z)continue;
      //const vec3 norm(det(defGradP[ps.p])*ps.norm.inner(inv(defGradP[ps.p])));
      //const vec3 b(curPosP[ps.p]);
      //const realT scale=.1*beamL; // length of gnuplot vectors
      //const vec3 e(scale*norm);
      //vf<<b<<tab<<e<<nwl;
   //}
   vf<<nwl<<nwl;
}

bool myPatch::afterStep(){
   history();
   //return false;
   return tmi->elapsedTime() < tmi->finalTime;
}

myPatch::~myPatch(){
   report.param("frameCount",frameCount);
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1){
      report.param("finalGap",gap0-abs(curPosP[i].y-refPosP[i].y));
   }
}

myPatch::myPatch(timeIntPtr&tm,
                 constitPtr&cst,
                 shapePtr&svar,
                 const realT bl,
                 const realT bh,
                 const realT bb):
                 tmi(tm),
                 shp(*svar),
                 beamL(bl),
                 beamH(bh),
                 beamB(bb){
   frameCount=0;
   ppd=comLineArg("ppd",vec3(2.,2.,2.));
   Ymod=comLineArg("Ymod",158.1e3);
   dens=comLineArg("dens",8.912e-3);
   pois=comLineArg("pois",.31);
   finalVolt=comLineArg("volt",100.);
   gap0=comLineArg("gap",.00375);
   const vec3 lo(0.,0.,0.);
   const vec3 hi(beamL,beamH,beamB);
   fillIndic(*this,*svar,lo,hi,boxIndic(lo,hi),vec3(2,2,2),false);
   cst=new updatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   //cst=new neoHookean(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("quasiCG")),0.,1.);
   for(indexT i=0;i<Npart();i+=1){
      massP[i]=refVolP[i]*dens;
   }
   partProbe(*svar,*this,vec3(.5*beamL*oneEps,0.,.5*beamB*oneEps));
   hf.open((string("history.xls")).c_str());  hf.precision(17);
   vf.open((string("vecPlot.xls")).c_str());  vf.precision(17);
   nf.open((string("nodes.xls")).c_str());  nf.precision(17);
   tf.open((string("tip.xls")).c_str());  tf.precision(17);
   history();
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   realT beamL=comLineArg("beamL",.4);
   realT beamH=comLineArg("beamH",.004);
   realT beamB=comLineArg("beamB",.12);
   int Nthick=comLineArg("Nthick",1);
   const realT cellThick=beamH/realT(Nthick);
   int Nlong=int(round(beamL/cellThick));
   int Nwide=int(round((beamB+2.*beamH)/cellThick));
   const vec3 lo(0.,-beamH,-beamH);
   const vec3 hi(beamL,beamH+beamH,beamB+beamH);
   makeShape(shp,Nlong,3*Nthick,Nwide,lo,hi,comLineArg("shape",string("MPM")));
   pch=new myPatch(tmi,cst,shp,beamL,beamH,beamB);
}










