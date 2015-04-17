// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

void makeFit(const vector<realT>&T,const vector<realT>&V,const string s); // defined in fit.cpp

struct exact;

struct myPatch:public patch{
   vector<realT>serTime,serVel;
   const shapeSC&shp;
   const realT beamL,beamH,beamB,Ymod,dens,pois,load;
   exact*ex;
   updatedElastic*cstPtr;
   ofstream hf,tf;
   const int mode;
   int frameCount,probeIdx;
   void vtkascii();
   bool afterStep();
   void gridVelocityBC(nodeArray<vec3>&,const realT&)const;
   myPatch(constitPtr&,shapeSC&,const realT,const realT,const realT);
   ~myPatch();
};

struct exact{
   const myPatch&cs;
   const realT L,H,B;
   exact(const myPatch&p,const realT l,const realT h,const realT b):cs(p),L(l),H(h),B(b){}
   realT kappa(const int i){
      const realT kap[]={4.7300407448627040,7.8532046240958376,10.995607838001671}; // fixed-fixed
      //const realT kap[]={1.875104069,4.694091133,7.854757438,10.99554073,14.13716839,17.27875953}; // cantilever
      if(i<1||i>3)throw stop("natFreq:bad natural frequency requested");
      return kap[i-1];
   }
   realT natFreq(const int i){
      const realT k=kappa(i);
      const realT Yeff=cs.Ymod/(1.-cs.pois*cs.pois); // plate correction
      return k*k*H/(L*L)*sqrt(Yeff/(12.*cs.dens));
   }
   realT modeShape(const int i,const realT&xBench){ // fixed-fixed
      const realT k=kappa(i);
      const realT fac=(sin(k)+sinh(k))/(cos(k)-cosh(k));
      const realT x=xBench/L;
      //return.5*(cos(k*x)-cosh(k*x)+fac*(sin(k*x)-sinh(k*x)));
      return cos(k*x)-cosh(k*x)+fac*(sin(k*x)-sinh(k*x));
   }
   //realT modeShape(const int i,const realT&xBench){ // cantilever
      //const realT k=kappa(i);
      //const realT fac=(cos(k)+cosh(k))/(sin(k)+sinh(k));
      //const realT x=xBench/L;
      //return.5*(cos(k*x)-cosh(k*x)-fac*(sin(k*x)-sinh(k*x)));
   //}
};

void myPatch::gridVelocityBC(nodeArray<vec3>&v,const realT&)const{
   for(int i=0;i<Nnode();++i){
      if(shp.gridPosN(i).x<machEps*beamL || shp.gridPosN(i).x>oneEps*beamL)v[i]=0.;
   }
}

void myPatch::vtkascii(){
   const string vtkName=string("zPart")+toStr(frameCount)+string(".vtk");
   ofstream vtkFile(vtkName.c_str());
   vtkFile<<"# vtk DataFile Version 3.0"<<nwl;
   vtkFile<<"Philip Wallstedt 3D MPM/GIMP"<<nwl;
   vtkFile<<"ASCII"<<nwl;
   vtkFile<<"DATASET POLYDATA"<<nwl;
   vtkFile<<"POINTS "<<Npart()<<" double"<<nwl;
   for(int i=0;i<Npart();++i)vtkFile<<curPosP[i]<<nwl;
   vtkFile<<nwl;
   vtkFile<<"POINT_DATA "<<Npart()<<nwl;
   vtkFile<<nwl;
   vtkFile<<"SCALARS vonMises double"<<nwl;
   vtkFile<<"LOOKUP_TABLE default"<<nwl;
   for(int i=0;i<Npart();++i)vtkFile<<vonMises(cstPtr->stress[i])<<nwl;
   vtkFile<<nwl;
   vtkFile<<"VECTORS displacement double"<<nwl;
   for(int i=0;i<Npart();++i)vtkFile<<curPosP[i]-refPosP[i]<<nwl;
}

bool myPatch::afterStep(){
   cerr.precision(14);
   //cerr<<"-- time   "<<elapsedTime<<endl;
   //cerr<<"refPos    "<<refPosP[probeIdx]<<endl;
   //cerr<<"xMP       "<<curPosP[probeIdx]<<endl;
   //cerr<<"velMP     "<<velP[probeIdx]<<endl;
   //cerr<<"massMP    "<<massP[probeIdx]<<endl;
   //cerr<<"dispMP    "<<curPosP[probeIdx]-refPosP[probeIdx]<<endl;
   //cerr<<"stressMP  "<<stressP[probeIdx]<<endl;
   tf<<elapsedTime<<tab<<velP[probeIdx].y<<tab<<load*cos(ex->natFreq(mode)*elapsedTime)*ex->modeShape(mode,curPosP[probeIdx].x)<<nwl;
   serTime.push_back(elapsedTime);
   serVel.push_back(velP[probeIdx].y);
   if(incCount%100==0){
      //vtkascii();
      hf<<'#'<<frameCount<<'\n';
      //for(int i=0;i<Nnode();++i){hf<<curPosN[i].x<<tab<<velIncN[i].y<<tab<<0.<<nwl;}
      hf<<'\n'<<endl;
      ++frameCount;
   }
   return elapsedTime<finalTime;
   //return incCount<1;
   //return false;
}

myPatch::~myPatch(){
   realT li=0.,l1=0.,l2=0.;
   int ct=0;
   for(int i=0;i<Npart();++i){
      ++ct;
      const realT uExact=load*cos(ex->natFreq(mode)*elapsedTime)*ex->modeShape(mode,curPosP[probeIdx].x);
      const realT uActual=velP[i].y;
      const realT er=abs(uActual-uExact);
      if(er>li)li=er;
      l1+=er;
      l2+=er*er;
   }
   l1=l1/realT(ct);
   l2=sqrt(l2/realT(ct));
   report.param("frameCount",frameCount);
   report.param("L1norm",l1);
   report.param("L2norm",l2);
   report.param("Linfnorm",li);
   makeFit(serTime,serVel,comLineArg("testName",string("run"))+string(".fit"));
}

myPatch::myPatch(constitPtr&cst,
                 shapeSC&svar,
                 const realT bl,
                 const realT bh,
                 const realT bb):
                 patch(svar),
                 shp(svar),
                 beamL(bl),
                 beamH(bh),
                 beamB(bb),
                 Ymod(comLineArg("Ymod",158.1e3)),
                 dens(comLineArg("dens",8.912e-3)),
                 pois(comLineArg("pois",.31)),
                 load(comLineArg("load",1.)),
                 mode(comLineArg("mode",1)){
   frameCount=probeIdx=0;
   cst=cstPtr=new updatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   ex=new exact(*this,beamL,beamH,beamB); // requires final Ymod and Dens
   elapsedTime=0.;
   report.param("EulerFreq",ex->natFreq(mode));
   finalTime=.00008;
   report.param("finalTime",finalTime);
   const vec3 lo(0.,0.,0.);
   const vec3 hi(beamL,beamH,beamB);
   boxIndic indic(lo,hi);
   fillIndic(*this,svar,lo,hi,indic,vec3(1,1,1));
   for(int i=0;i<Npart();i+=1){
      velP[i]=vec3(0.,load*ex->modeShape(mode,refPosP[i].x),0.);
      massP[i]=refVolP[i]*dens;
   }
   probeIdx=partProbe(*this,vec3(.01,.002,.002));
   hf.open("history.xls");  hf.precision(17);
   tf.open("tip.xls");      tf.precision(17);
   tf<<elapsedTime<<tab<<velP[probeIdx].y<<nwl;
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   const realT beamL=comLineArg("beamL",.02);
   const realT beamH=comLineArg("beamH",.004);
   const realT beamB=comLineArg("beamB",.004);
   const vec3 lo(0,0,0);
   const vec3 hi(beamL,beamH,beamB);
   makeShape(shp,5,1,1,lo,hi,comLineArg("shape",string("MPM")));
   pch=new myPatch(cst,*shp,beamL,beamH,beamB);
   makeTimeInt(tmi,pch,cst,shp,comLineArg("tInt",string("CD")));
   cerr<<endl;
}












