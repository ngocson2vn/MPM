// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

class hollowBox{
   const vec3 orig,half;
   boxIndic outside,inside;
public:
   hollowBox(const vec3&l,
             const vec3&h,
             const realT&thick):
                orig((h+l)*.5),
                half((h-l)*.5),
                outside(orig-half,orig+half),
                inside(orig-half+thick,orig+half-thick,true){}
   bool operator()(const vec3&pt)const{return outside(pt) && inside(pt);}
};

struct myPatch:public patch{
   mutable partArray<vec3>surfaceP;
   const realT Ymod,dens,pois,pressure;
   indexT frameCount;
   updatedElastic*nh;
   timeIntPtr&tmi;
   shapePtr&shp;
   void vtkascii();
   bool afterStep();
   void gridVelocityBC(nodeArray<vec3>&,const realT&)const;
   void makeExternalForces(nodeArray<vec3>&,const realT&)const;
   myPatch(timeIntPtr&tmi,constitPtr&cst,shapePtr&s);
   ~myPatch();
};

void myPatch::vtkascii(){
   const string vtkName=string("zPart")+toStr(tmi->incCount())+string(".vtk");
   ofstream vtkFile(vtkName.c_str());
   vtkFile<<"# vtk DataFile Version 3.0"<<nwl;
   vtkFile<<"Philip Wallstedt MPM-GIMP"<<nwl;
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

void myPatch::makeExternalForces(nodeArray<vec3>&gravityN,const realT&delt)const{
   const realT loadFactor=(tmi->elapsedTime() - tmi->timeStep() + delt) / tmi->finalTime;
   report.param("loadFactor",loadFactor);
   const vec3 snapNorm(-1.,-1.,-1.);
   for(indexT i=0;i<Npart();++i)if(surfNormP[i].inner(snapNorm) > .707){
      const realT refArea=surfAreaP[i].inner(surfNormP[i]);
      surfaceP[i]=loadFactor*pressure*det(defGradP[i])*refArea*surfNormP[i].inner(inv(defGradP[i]));
   }
   shp->integrate(surfaceP,gravityN);
   for(indexT i=0;i<Nnode();++i)gravityN[i]/=massN[i];
}

void myPatch::gridVelocityBC(nodeArray<vec3>&v,const realT&)const{
   for(indexT i=0;i<Nnode();++i){
      if(shp->gridPosN(i).x<machEps)v[i].x=0.;
      if(shp->gridPosN(i).y<machEps)v[i].y=0.;
      if(shp->gridPosN(i).z<machEps)v[i].z=0.;
   }
}

bool myPatch::afterStep(){
   vtkascii();
   cerr<<'.';
   return tmi->elapsedTime() < tmi->finalTime;
   //return false;
}

myPatch::~myPatch(){
   indexT count=1;
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1){
      string tag=string("probe ")+toStr(count)+string(" on rank ")+toStr(cartRank);
      report.param(tag.c_str(),mag(curPosP[i]-refPosP[i]));
      cerr<<tag<<tab<<mag(curPosP[i]-refPosP[i])<<endl;
      ++count;
   }
}

myPatch::myPatch(timeIntPtr&tm,
                 constitPtr&cst,
                 shapePtr&svar):
                    Ymod(comLineArg("Ymod",210e9)),
                    dens(comLineArg("dens",7850.)),
                    pois(comLineArg("pois",.31)),
                    pressure(comLineArg("pres",1.)),
                    tmi(tm),
                    shp(svar){
   const vec3 domLo(0.,0.,0.);
   const vec3 domHi(5.,5.,5.);
   int Nc=comLineArg("Ncell",5);
   makeShape(shp,Nc,Nc,Nc,domLo,domHi,comLineArg("shape",string("MPM")));
   const vec3 boxLo(-4.,-4.,-4.); // only using first quadrant
   const vec3 boxHi( 4., 4., 4.);
   fillIndic(*this,*svar,boxLo,boxHi,hollowBox(boxLo,boxHi,1.),vec3(2,2,2),false);
   cst=nh=new updatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("quasiCG")),0.,1.);
   for(indexT i=0;i<Npart();i+=1){
      massP[i]=refVolP[i]*dens;
   }
   partProbe(*svar,*this,oneEps*vec3(4.,4.,4.));
   partProbe(*svar,*this,oneEps*vec3(4.,4.,0.));
   partProbe(*svar,*this,oneEps*vec3(0.,0.,4.));
   partProbe(*svar,*this,oneEps*vec3(2.,2.,4.));
   vtkascii();
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){pch=new myPatch(tmi,cst,shp);}


