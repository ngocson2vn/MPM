// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef VIRTUAL_INTERFACE_H
#define VIRTUAL_INTERFACE_H

#include<cstdlib>
#include<typeinfo>
#include "shapePre.h"
//using namespace std;

struct patch{
   partArray<mat3>velGradP,defGradP,stressP;
   partArray<vec3>refPosP,curPosP,velP,halfLenP,refHalfP,surfNormP,surfAreaP;
   partArray<realT>massP,volumeP,refVolP;
   partArray<indexT>probeP;
   nodeArray<vec3>velN,fintN,bodyAccelN,velIncN;
   nodeArray<realT>massN;
   patch(){for(indexT i=0;i<Npart();++i)probeP[i]=0;}
   virtual ~patch(){}
   virtual bool afterStep()=0;
   virtual void gridVelocityBC(nodeArray<vec3>&,const realT&)const{};
   virtual void makeExternalForces(nodeArray<vec3>&,const realT&)const{};
};
typedef patch*patchPtr;

// Every type of shape function must provide this
struct shapeSC:public shapePre{
   shapeSC(const indexT Nx,const indexT Ny,const indexT Nz,const vec3 lo,const vec3 hi):shapePre(Nx,Ny,Nz,lo,hi){}
   void integrate(const partArray<realT>&pu,nodeArray<realT>&gu)const;
   void integrate(const partArray<vec3>&pu,nodeArray<vec3>&gu)const;
   void divergence(const partArray<mat3>&stressP,const partArray<realT>&volumeP,nodeArray<vec3>&forceN)const;
   void interpolate(partArray<vec3>&pu,const nodeArray<vec3>&gu)const;
   void gradient(partArray<mat3>&pu,const nodeArray<vec3>&gu)const;
   void gradient(partArray<vec3>&pu,const nodeArray<realT>&gu)const;
   virtual void updateContribList(const patch&)=0;
   virtual~shapeSC(){}
};
typedef shapeSC*shapePtr;
void makeShape(shapePtr&shp,const indexT Nx,const indexT Ny,const indexT Nz,const vec3 lo,const vec3 hi,string s);

// constitutive super class
struct constitutiveSC{
   virtual void update(const realT&)=0;
   __attribute__((noreturn)) virtual void revert(){throw stop("constitutiveSC: revert() and save() must be defined if using implicit methods");}
   __attribute__((noreturn)) virtual void save()  {throw stop("constitutiveSC: revert() and save() must be defined if using implicit methods");}
   virtual realT waveSpeed()const=0;
   virtual ~constitutiveSC(){}
};
typedef constitutiveSC*constitPtr;

// Several algorithmic variations are used, but they all provide these functions
class timeIntSC{
   realT _timeStep,nominalStep,Nstep,_elapsedTime;
   indexT _incCount;
public:
   const realT initialTime,finalTime;
   realT elapsedTime()const{return _elapsedTime;}
   indexT incCount(){return _incCount;}
   timeIntSC(const patch&pch,const constitutiveSC&cst,const shapeSC&shp,realT&t0,realT&tf);
   virtual~timeIntSC(){}
   virtual void advance(const shapeSC&shp)=0;
   void nextStep(const patch&pch);
   realT timeStep()const{return _timeStep;}
   indexT intEstSteps()const{return indexT(round(Nstep));}
};
typedef timeIntSC*timeIntPtr;
void makeTimeInt(timeIntPtr&ti,patch*pch,constitutiveSC*cst,shapeSC*shp,string s,realT t0,realT tf=huge);

// io.cpp - utility functions for setting up arrangements of particles
bool partProbe(const shapeSC&shp,patch&pch,const vec3&pos);
indexT nodeProbe(const shapeSC&shp,const vec3&pos);

// Every input file must define initRun
void initRun(patchPtr&,constitPtr&,timeIntPtr&,shapePtr&);







#endif
