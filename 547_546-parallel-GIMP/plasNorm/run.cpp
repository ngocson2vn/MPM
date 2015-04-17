// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

struct myPatch:public patch{
   mutable partArray<vec3>bodyP;
   ofstream hf;
   realT load;
   int frameCount;
   timeIntPtr&tmi;
   const shapeSC&shp;
   J2plasticLinearIsoKin*nh;
   bool afterStep();
   void vtkascii();
   void history();
   void gridVelocityBC(nodeArray<vec3>&,const realT&)const;
   void makeExternalForces(nodeArray<vec3>&,const realT&)const;
   myPatch(timeIntPtr&,constitPtr&cst,shapePtr&s);
   ~myPatch();
};

void myPatch::makeExternalForces(nodeArray<vec3>&gravityN,const realT&)const{
   for(indexT i=0;i<Npart();++i)bodyP[i]=vec3(massP[i] * tmi->elapsedTime()*load,0.,0.);
   shp.integrate(bodyP,gravityN);
   for(indexT i=0;i<Nnode();++i)gravityN[i]/=massN[i];
}

void myPatch::gridVelocityBC(nodeArray<vec3>&v,const realT&)const{
   for(indexT i=0;i<Nnode();++i){
      if(shp.gridPosN(i).x<machEps)v[i].x=0.;
   }
}

void myPatch::vtkascii(){
   const string vtkName=string("zPart")+toStr(tmi->incCount())+string(".vtk");
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
   for(indexT i=0;i<Npart();++i)vtkFile<<vonMises(nh->stress[i])<<nwl;
   vtkFile<<nwl;
   vtkFile<<"VECTORS displacement double"<<nwl;
   for(indexT i=0;i<Npart();++i)vtkFile<<curPosP[i]-refPosP[i]<<nwl;
}

void myPatch::history(){
   frameCount+=1;
   hf<<'#'<<frameCount<<'\n';
   for(indexT i=0;i<Npart();++i){
      if(refPosP[i].z>.5&&refPosP[i].z<.5+.5*shp.cell.z){
         if(refPosP[i].y>.5)hf<<curPosP[i]<<tab<<vonMises(nh->stress[i])<<nwl;
         else               hf<<curPosP[i]<<tab<<trace(nh->stress[i])/3.<<nwl;
      }
   }
   hf<<"\n\n";
}

bool myPatch::afterStep(){
   vtkascii();
   if(tmi->incCount()%1==0)history();
   return tmi->elapsedTime() < tmi->finalTime;
   //return false;
}

myPatch::~myPatch(){
   report.param("frameCount",frameCount);
   hf.close();
}

myPatch::myPatch(timeIntPtr&tm,constitPtr&cst,shapePtr&svar):tmi(tm),shp(*svar){
   frameCount=0;
   load=comLineArg("load",1.);    report.param("load",load);
   const vec3 lo(0,0,0);
   const vec3 hi(1,1,1);
   fillIndic(*this,*svar,lo,hi,boxIndic(lo,hi),vec3(2,2,2),false);
   const realT Ymod=comLineArg("Ymod",200e9);
   const realT pois=comLineArg("pois",.3);
   const realT dens=comLineArg("dens",8100.);
   //cst=nh=new J2plasticLinearIsoKin(Ymod,pois,dens,1e9,11e9,0.  ,velGradP,defGradP,volumeP,stressP); // iso
   cst=nh=new J2plasticLinearIsoKin(Ymod,pois,dens,1e9,0.  ,11e9,velGradP,defGradP,volumeP,stressP); // kin
   //cst=nh=new J2plasticLinearIsoKin(Ymod,pois,dens,1e9,0.  ,0.  ,velGradP,defGradP,volumeP,stressP); // perf
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("CD")),0.,.001);
   for(indexT i=0;i<Npart();i+=1){
      massP[i]=refVolP[i]*dens;
   }
   hf.open("history.xls");  hf.precision(14);
   history();
   vtkascii();
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   const vec3 lo(0.,-1.,-1.);
   const vec3 hi(2.,2.,2.);
   int Nc=comLineArg("Ncell",10);
   makeShape(shp,2*Nc,3*Nc,3*Nc,lo,hi,comLineArg("shape",string("GIMP")));
   pch=new myPatch(tmi,cst,shp);
}


