// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

struct myPatch:public patch{
   mutable partArray<vec3>surfaceP;
   vector<realT>serTime,serVel;
   timeIntPtr&tmi;
   const shapeSC&shp;
   const realT beamL,beamH,beamB,Ymod,dens,pois,load;
   updatedElastic*cstPtr;
   ofstream hf,tf;
   int frameCount;
   realT EBdef(const realT&x){ // fixed-fixed beam with uniform load
      const realT momI=beamB*beamH*beamH*beamH/12.;
      const realT w=-load*beamB; // load per length
      const realT Yeff=Ymod/(1.-pois*pois);
      return w*x*x*(beamL-x)*(beamL-x)/(24.*Yeff*momI);
   }
   //realT defGivenPres(const realT fperA){                    // fixed-fixed beam with uniform load
      //const realT I=beamB*beamH*beamH*beamH/12.;
      //const realT w=-fperA*beamB;                            // load per length
      //const realT x=.5*beamL;                                // midpoint of beam
      //return w*x*x*(beamL-x)*(beamL-x)/(24.*Ymod*I);
   //}
   //realT presGivenFracDef(const realT fracDef){              // fixed-fixed beam with uniform load
      //const realT I=beamB*beamH*beamH*beamH/12.;
      //const realT deflect=fracDef*beamH/beamB;               // find load required to deflect by this amount
      //const realT x=.5*beamL;                                // midpoint of beam
      //return 24.*Ymod*I*deflect/(x*x*(beamL-x)*(beamL-x));
   //}
   void vtkascii();
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

void myPatch::makeExternalForces(nodeArray<vec3>&gravityN,const realT&delt)const{ // new surface
   const realT loadFactor=(tmi->elapsedTime()-tmi->timeStep()+delt)/tmi->finalTime;
   report.param("loadFactor",loadFactor);
   const vec3 snapNorm(0,1,0);
   for(indexT i=0;i<Npart();++i)if(surfNormP[i].y>0.){
      const realT currArea=surfAreaP[i].inner(snapNorm);
      surfaceP[i]=(-load*loadFactor*det(defGradP[i])*currArea*snapNorm.inner(inv(defGradP[i])));
   }
   shp.integrate(surfaceP,gravityN);
   for(indexT i=0;i<Nnode();++i)gravityN[i]/=massN[i];
}

//void myPatch::makeExternalForces(nodeArray<vec3>&gravityN,const realT&delt)const{ // body force
   //const realT loadFactor=(elapsedTime-timeStep+delt)/finalTime;
   //report.param("loadFactor",loadFactor);
   //const realT partForce=-loadFactor*load*beamL*beamB/realT(Npart());
   //for(indexT i=0;i<Npart();++i)surfaceP[i]=partForce;
   //shp.integrate(surfaceP,gravityN);
   //for(indexT i=0;i<Nnode();++i)gravityN[i]/=massN[i];
//}

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
   const realT loadFactor=tmi->elapsedTime()/tmi->finalTime;
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1){
      tf<<loadFactor<<tab<<curPosP[i].y-refPosP[i].y<<tab<<loadFactor*EBdef(refPosP[i].x)<<nwl;
      serTime.push_back(tmi->elapsedTime());
      serVel.push_back(curPosP[i].y-refPosP[i].y);
   }
   if(tmi->incCount()%1==0){
      //vtkascii();
      hf<<'#'<<frameCount<<'\n';
      //for(indexT i=0;i<Nnode();++i){hf<<curPosN[i].x<<tab<<velIncN[i].y<<tab<<0.<<nwl;}
      //for(indexT i=0;i<Npart();++i){hf<<curPosP[i]<<tab<<surfAreaP[i]<<nwl;}
      for(indexT i=0;i<Npart();++i){hf<<curPosP[i]<<tab<<surfAreaP[i]<<nwl;}
      hf<<'\n'<<endl;
      ++frameCount;
   }
   return tmi->elapsedTime() < tmi->finalTime;
   //return false;
}

myPatch::~myPatch(){
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1){
      report.param("frameCount",frameCount);
      report.param("computedDeflection",curPosP[i].y-refPosP[i].y);
      report.param("EulerBerDeflection",EBdef(refPosP[i].x));
      cerr<<"deflection "<<curPosP[i].y-refPosP[i].y<<endl;
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
                 beamB(bb),
                 Ymod(comLineArg("Ymod",158.1e3)),
                 dens(comLineArg("dens",8.912e-3)),
                 pois(comLineArg("pois",.31)),
                 load(comLineArg("load",1.)){
   frameCount=0;
   cst=cstPtr=new updatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("CD")),0.,1.);
   const vec3 lo(0.,0.,0.);
   const vec3 hi(beamL,beamH,beamB);
   fillIndic(*this,*svar,lo,hi,boxIndic(lo,hi),vec3(2,2,2),false);
   for(indexT i=0;i<Npart();i+=1){
      massP[i]=refVolP[i]*dens;
   }
   partProbe(*svar,*this,.5*vec3(beamL,beamH,beamB));
   hf.open("history.xls");  hf.precision(17);
   tf.open("tip.xls");      tf.precision(17);
   tf<<tmi->elapsedTime()<<tab<<0.<<tab<<0.<<nwl;
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   const realT beamL=comLineArg("beamL",.4);
   const realT beamH=comLineArg("beamH",.004);
   const realT beamB=comLineArg("beamB",.12);
   const int Nthick=comLineArg("Nthick",1);
   const realT cellThick=beamH/realT(Nthick);
   const int Nwide=comLineArg("Nwide",int(round(beamB/cellThick)));
   const realT cellWide=beamB/realT(Nwide);
   const int Nlong=int(round(beamL/cellThick));
   const int difWide=int(round(1.25*realT(Nwide)))-Nwide;
   const vec3 lo(0.   ,-2.*beamH,-cellWide*difWide);
   const vec3 hi(beamL, 3.*beamH, cellWide*(Nwide+difWide));
   makeShape(shp,Nlong,5*Nthick,Nwide+2*difWide,lo,hi,comLineArg("shape",string("MPM")));
   pch=new myPatch(tmi,cst,shp,beamL,beamH,beamB);
}












