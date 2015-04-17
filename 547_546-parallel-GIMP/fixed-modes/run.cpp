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
   timeIntPtr&tmi;
   shapePtr&shp;
   const realT beamL,beamH,beamB,Ymod,dens,pois,load;
   exact*ex;
   updatedElastic*cstPtr;
   ofstream hf,tf;
   const int mode;
   int frameCount;
   void vtkascii();
   bool afterStep();
   void gridVelocityBC(nodeArray<vec3>&,const realT&)const;
   myPatch(timeIntPtr&tmi,constitPtr&cst,shapePtr&s,const realT,const realT,const realT);
   ~myPatch();
};

struct exact{
   const myPatch&cs;
   const realT L,H,B;
   exact(const myPatch&p,const realT l,const realT h,const realT b):cs(p),L(l),H(h),B(b){}
   realT kappa(const indexT i){
      const realT kap[]={4.7300407448627040,7.8532046240958376,10.995607838001671}; // fixed-fixed
      //const realT kap[]={1.875104069,4.694091133,7.854757438,10.99554073,14.13716839,17.27875953}; // cantilever
      if(i<1||i>3)throw stop("natFreq:bad natural frequency requested");
      return kap[i-1];
   }
   realT natFreq(const indexT i){
      const realT k=kappa(i);
      const realT Yeff=cs.Ymod/(1.-cs.pois*cs.pois); // plate correction
      return k*k*H/(L*L)*sqrt(Yeff/(12.*cs.dens));
   }
   realT modeShape(const indexT i,const realT&xBench){ // fixed-fixed
      const realT k=kappa(i);
      const realT fac=(sin(k)+sinh(k))/(cos(k)-cosh(k));
      const realT x=xBench/L;
      return.5*(cos(k*x)-cosh(k*x)+fac*(sin(k*x)-sinh(k*x)));
   }
   //realT modeShape(const indexT i,const realT&xBench){ // cantilever
      //const realT k=kappa(i);
      //const realT fac=(cos(k)+cosh(k))/(sin(k)+sinh(k));
      //const realT x=xBench/L;
      //return.5*(cos(k*x)-cosh(k*x)-fac*(sin(k*x)-sinh(k*x)));
   //}
};

void myPatch::gridVelocityBC(nodeArray<vec3>&v,const realT&)const{
   for(indexT i=0;i<Nnode();++i){
      if(shp->gridPosN(i).x<machEps*beamL || shp->gridPosN(i).x>oneEps*beamL)v[i]=0.;
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
   for(indexT i=0;i<Npart();++i)vtkFile<<curPosP[i]<<nwl;
   vtkFile<<nwl;
   vtkFile<<"POINT_DATA "<<Npart()<<nwl;
   vtkFile<<nwl;
   vtkFile<<"SCALARS vonMises double"<<nwl;
   vtkFile<<"LOOKUP_TABLE default"<<nwl;
   for(indexT i=0;i<Npart();++i)vtkFile<<vonMises(cstPtr->stress[i])<<nwl;
   vtkFile<<nwl;
   vtkFile<<"VECTORS displacement double"<<nwl;
   for(indexT i=0;i<Npart();++i)vtkFile<<curPosP[i]-refPosP[i]<<nwl;
}

bool myPatch::afterStep(){
   if(tmi->incCount()%1000==0)for(indexT i=0;i<Npart();++i)if(probeP[i]==1)cerr<<cartRank<<" probe "<<refPosP[i]<<"   disp"<<tab<<curPosP[i]-refPosP[i]<<endl;
   serTime.push_back(tmi->elapsedTime());
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1){
      tf<<tmi->elapsedTime()<<tab<<velP[i].y<<tab<<load*cos(ex->natFreq(mode)*tmi->elapsedTime())*ex->modeShape(mode,curPosP[i].x)<<nwl;
      serVel.push_back(velP[i].y);
   }
   if(tmi->incCount()%100==0){
      //vtkascii();
      hf<<'#'<<frameCount<<'\n';
      //for(indexT i=0;i<Nnode();++i){hf<<curPosN[i].x<<tab<<velIncN[i].y<<tab<<0.<<nwl;}
      hf<<'\n'<<endl;
      ++frameCount;
   }
   return tmi->elapsedTime() < tmi->finalTime;
   //return false;
}

myPatch::~myPatch(){
   realT li=0.,l1=0.,l2=0.;
   int ct=0;
   realT uExact=0.;
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1)uExact=load*cos(ex->natFreq(mode)*tmi->elapsedTime())*ex->modeShape(mode,curPosP[i].x);
   for(indexT i=0;i<Npart();++i){
      ++ct;
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

myPatch::myPatch(timeIntPtr&tm,
                 constitPtr&cst,
                 shapePtr&svar,
                 const realT bl,
                 const realT bh,
                 const realT bb):
                 tmi(tm),
                 shp(svar),
                 beamL(bl),
                 beamH(bh),
                 beamB(bb),
                 Ymod(comLineArg("Ymod",158.1e3)),
                 dens(comLineArg("dens",8.912e-3)),
                 pois(comLineArg("pois",.31)),
                 load(comLineArg("load",1.)),
                 mode(comLineArg("mode",1)){


   const int Nthick=comLineArg("Nthick",1);
   const realT cellThick=beamH/realT(Nthick);
   const int Nwide=comLineArg("Nwide",int(round(beamB/cellThick)));
   const realT cellWide=beamB/realT(Nwide);
   const int Nlong=int(round(beamL/cellThick));
   const int difWide=int(round(1.25*realT(Nwide)))-Nwide;
   //cerr<<"Nthick  "<<Nthick<<endl;
   //cerr<<"Nwide   "<<Nwide<<endl;
   //cerr<<"Nlong   "<<Nlong<<endl;
   //cerr<<"difWide "<<difWide<<endl;
   const vec3 dlo(0.   ,-2.*beamH,-cellWide*difWide);
   const vec3 dhi(beamL, 2.*beamH, cellWide*(Nwide+difWide));
   makeShape(shp,Nlong,5*Nthick,Nwide+2*difWide,dlo,dhi,comLineArg("shape",string("GIMP")));

   frameCount=0;
   cst=cstPtr=new updatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   ex=new exact(*this,beamL,beamH,beamB); // requires final Ymod and Dens
   report.param("EulerFreq",ex->natFreq(mode));
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("CD")),0.,4./(ex->natFreq(mode)/(2.*pi)));
   report.param("finalTime",tmi->finalTime);
   const vec3 lo(0.,0.,0.);
   const vec3 hi(beamL,beamH,beamB);
   boxIndic indic(lo,hi);
   fillIndic(*this,*svar,lo,hi,indic);
   for(indexT i=0;i<Npart();i+=1){
      velP[i]=vec3(0.,load*ex->modeShape(mode,refPosP[i].x),0.);
      massP[i]=refVolP[i]*dens;
   }
   partProbe(*shp,*this,.5*vec3(beamL,beamH,beamB));
   hf.open("history.xls");  hf.precision(17);
   tf.open("tip.xls");      tf.precision(17);
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1)tf<<tmi->elapsedTime()<<tab<<velP[i].y<<nwl;
   cerr<<endl;
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   const realT beamL=comLineArg("beamL",.4);
   const realT beamH=comLineArg("beamH",.004);
   const realT beamB=comLineArg("beamB",.12);
   pch=new myPatch(tmi,cst,shp,beamL,beamH,beamB);
}












