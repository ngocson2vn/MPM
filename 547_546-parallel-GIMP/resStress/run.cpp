// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

struct myPatch:public patch{
   mutable partArray<vec3>forceP;
   mutable partArray<mat3>initStsP;
   mutable nodeArray<vec3>finitN;
   timeIntPtr&tmi;
   const shapeSC&shp;
   const realT beamL,beamH,beamB,Ymod,dens,pois,force0,forceM,resStress;
   updatedElastic*cstPtr;
   ofstream hf,tf;
   int frameCount;
   indexT padSize0,padSizeM;
   realT defPoint(const realT x,const realT xforce,const realT force){
      const realT momI=beamB*beamH*beamH*beamH/12.;
      const realT a=xforce;
      const realT b=beamL-xforce;
      const realT t1=force/(6.*(Ymod/(1.-pois*pois))*momI);
      //const realT t1=force/(6.*Ymod*momI);
      const realT t2=b*b*x*x*x*(beamL+2.*a)/(beamL*beamL*beamL);
      const realT t3=3.*a*b*b*x*x/(beamL*beamL);
      const realT xa=x-a;
      const realT t4=(x>a?xa*xa*xa:0.);
      return t1*(t2-t3-t4);
   }
   realT defPad(const indexT&id){
      realT def=0.;
      if(id==1)for(indexT i=0;i<Npart();++i)if(probeP[i]==id)def+=(curPosP[i].y-refPosP[i].y)/realT(padSize0);
      if(id==2)for(indexT i=0;i<Npart();++i)if(probeP[i]==id)def+=(curPosP[i].y-refPosP[i].y)/realT(padSizeM);
      MPI_Allreduce(MPI_IN_PLACE,&def,1,MPI_DOUBLE,MPI_SUM,cartComm);
      return def;
   }
   bool afterStep();
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
   for(indexT i=0;i<Npart();++i)forceP[i]=0.;
   for(indexT i=0;i<Npart();++i){
      if(probeP[i]==1)forceP[i].y+=force0/realT(padSize0);
      if(probeP[i]==2)forceP[i].y+=forceM/realT(padSizeM);
   }
   shp.integrate(forceP,gravityN);

   //const realT loadFactor=(incCount==0?0.:(elapsedTime-delt)/(finalTime-delt)); // zero for first step, one for last step
   //const realT loadFactor=(elapsedTime+delt)/finalTime;

   //realT loadFactor=elapsedTime/(finalTime-3.*delt); // equilibrate for last three steps
   //if(elapsedTime==0.)loadFactor=0.;
   //if((elapsedTime+machEps*finalTime)>finalTime-3.*delt)loadFactor=1.;
   //if(abs(elapsedTime+delt-finalTime)<machEps)loadFactor=1.;

   const realT loadFactor=(tmi->incCount()>0?1.:0.); // apply force, apply stress, equilibrate for remaining steps

   cerr<<'L'<<loadFactor<<' ';
   report.param("loadFactor",loadFactor);
   for(indexT i=0;i<Npart();++i)initStsP[i]=mat3(loadFactor*resStress,0,0,0,0,0,0,0,0);
   shp.divergence(initStsP,volumeP,finitN);
   for(indexT i=0;i<Nnode();++i)gravityN[i]+=finitN[i];
   for(indexT i=0;i<Nnode();++i)gravityN[i]/=massN[i];
}

bool myPatch::afterStep(){
   //const realT loadFactor=tmi->elapsedTime()/tmi->finalTime;
   //tf<<loadFactor<<tab
     //<<defPad(2)<<tab
     //<<defPad(1)<<tab
     //<<-defPoint(.5*beamL,.5*beamL,forceM*loadFactor)<<nwl;
   //cerr<<"dM d0 "<<defPad(2)<<tab<<defPad(1)<<endl;
   hf<<'#'<<frameCount<<'\n';
   for(indexT i=0;i<Npart();++i)hf<<curPosP[i]<<tab<<velP[i]<<nwl;
   hf<<nwl<<endl;
   ++frameCount;
   return tmi->elapsedTime() < tmi->finalTime;
   //return false;
}

myPatch::~myPatch(){
   report.param("frameCount",frameCount);
   report.param("exact midpoint deflection",defPoint(.5*beamL,.5*beamL,forceM));
   report.param("exact midpoint stiffness",forceM/defPoint(.5*beamL,.5*beamL,forceM));
   report.param("final stiffness pad0",force0/defPad(1));
   report.param("final stiffness padM",forceM/defPad(2));
   report.param("deflectPad0",defPad(1));
   report.param("deflectPadM",defPad(2));
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
                 beamB(bb),
                 Ymod(comLineArg("Ymod",158.1e3)),
                 dens(comLineArg("dens",8.912e-3)),
                 pois(comLineArg("pois",.31)),
                 force0(comLineArg("f0",0.)),
                 forceM(comLineArg("fM",0.)),
                 resStress(comLineArg("RS",0.)){
   frameCount=0;
   padSize0=padSizeM=0;
   cst=cstPtr=new updatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("quasiCG")),0.,1.);
   const vec3 lo(0.,0.,0.);
   const vec3 hi(beamL,beamH,beamB);
   fillIndic(*this,*svar,lo,hi,boxIndic(lo,hi),vec3(2,2,2),false);
   for(indexT i=0;i<Npart();i+=1){
      massP[i]=refVolP[i]*dens;
   }
   const vec3 orig0(.5*beamL,beamH,shp.cell.z); // near edge
   const vec3 origM(.5*beamL,beamH,.5*beamB     ); // center
   for(indexT i=0;i<Npart();++i)if(surfNormP[i].y>=.999){
      const vec3&p=refPosP[i];
      if(abs(p.x-orig0.x)<.5*shp.cell.x && abs(p.z-orig0.z)<.5*shp.cell.z){cerr<<"pad0 "; probeP[i]=1; ++padSize0;}
      if(abs(p.x-origM.x)<.5*shp.cell.x && abs(p.z-origM.z)<.5*shp.cell.z){cerr<<"padM "; probeP[i]=2; ++padSizeM;}
   }
   cerr<<endl;
   MPI_Allreduce(MPI_IN_PLACE,&padSize0,1,MPI_UNSIGNED,MPI_SUM,cartComm);
   MPI_Allreduce(MPI_IN_PLACE,&padSizeM,1,MPI_UNSIGNED,MPI_SUM,cartComm);
   cerr<<"pad sizes "<<padSize0<<tab<<padSizeM<<endl;
   hf.open("history.xls");  hf.precision(17);
   tf.open("tip.xls");      tf.precision(17);
   tf<<tmi->elapsedTime()<<tab<<0.<<tab<<0.<<tab<<0.<<nwl;
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   const realT beamL=comLineArg("beamL",.4);
   const realT beamH=comLineArg("beamH",.004);
   const realT beamB=comLineArg("beamB",.12);
   const int Nthick=comLineArg("Nthick",1);
   const realT cellThick=beamH/realT(Nthick);
   const int Nwide=comLineArg("Nwide",int(round(beamB/cellThick)));
   const int Nlong=int(round(beamL/cellThick));
   const vec3 lo;
   const vec3 hi(beamL,beamH,beamB);
   makeShape(shp,Nlong,Nthick,Nwide,lo,hi,comLineArg("shape",string("MPM")));
   pch=new myPatch(tmi,cst,shp,beamL,beamH,beamB);
}












