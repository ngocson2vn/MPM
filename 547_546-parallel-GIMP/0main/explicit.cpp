// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include "timeInt.h"
using namespace std;

//void timeIntOps::setMinMass(const patch&pch){
   //realT sum=0.;
   //for(indexT i=0;i<Npart();++i)sum+=pch.massP[i];
   //MPI_Allreduce(MPI_IN_PLACE,&sum,1,MPI_DOUBLE,MPI_SUM,cartComm);
   //minMass=machEps*sum;
   //report.param("minMass",minMass);
//}

void timeIntOps::setMinMass(const patch&pch){
   const realT*aptr=&pch.massP[0];
   complex<realT>sum(0.,0.);
   complex<realT>addend(0.,0.);
   for(indexT j=0;j<Npart();++j){
      addend.real()=aptr[j];
      DDadd(&addend.real(),&sum.real(),NULL,NULL);
   }
   MPI_Allreduce(MPI_IN_PLACE,&sum,1,MPI::DOUBLE_COMPLEX,MPI_DDadd,cartComm);
   minMass=machEps*sum.real();
   report.param("minMass",minMass);
}

void timeIntOps::integrateMass(const partArray<realT>&massP,nodeArray<realT>&massN){
   shp.integrate(massP,massN);
   for(indexT i=0;i<Nnode();++i)if(massN[i]<minMass)massN[i]=minMass;
}

void timeIntOps::solveMomentum(const patch&pch,
                               const realT&dt,
                               nodeArray<vec3>&fintN,
                               nodeArray<vec3>&velN,
                               nodeArray<vec3>&velIncN){
   shp.divergence(pch.stressP,pch.volumeP,fintN);
   for(indexT i=0;i<Nnode();++i)velIncN[i]=(-velN[i]);
   for(indexT i=0;i<Nnode();++i)velN[i]+=(pch.bodyAccelN[i]+fintN[i]/pch.massN[i])*dt;
   pch.gridVelocityBC(velN,dt);
   for(indexT i=0;i<Nnode();++i)velIncN[i]+=velN[i];
}

void timeIntOps::advancePosition(const nodeArray<vec3>&velN,const realT&dt,partArray<vec3>&curPosP){
   shp.interpolate(incP,velN);
   for(indexT i=0;i<Npart();++i)curPosP[i]+=incP[i]*dt;
}

void timeIntOps::advanceVelocity(const nodeArray<vec3>&velIncN,const realT&,partArray<vec3>&velP){
   shp.interpolate(incP,velIncN);
   for(indexT i=0;i<Npart();++i)velP[i]+=incP[i];
}

void timeIntOps::mapVelPtoVelN(const patch&pch,const realT&dt,const partArray<vec3>&velP,nodeArray<vec3>&velN){
   for(indexT i=0;i<Npart();++i)momP[i]=velP[i]*pch.massP[i];
   shp.integrate(momP,momN);
   for(indexT i=0;i<Nnode();++i)velN[i]=momN[i]/pch.massN[i];
   pch.gridVelocityBC(velN,dt);
}

void timeIntOps::advanceConstit(patch&pch,constitutiveSC&cst,const realT&dt){
   shp.gradient(pch.velGradP,pch.velN);
   cst.update(dt);
   for(indexT i=0;i<Npart();++i)pch.halfLenP[i]=pch.refHalfP[i]*diag(pch.defGradP[i]);
}

void makeTimeInt(timeIntPtr&ti,patch*pch,constitutiveSC*cst,shapeSC*shp,string s,realT t0,realT tf){
   if(cst==NULL)throw stop("constitutiveSC pointer cst must be set in makePatch");
   const bool stat=true;
   const bool sym=true;
   const bool remap=true;
   if(s=="momentum"){ti=new explicitUSL< remap>                        (*pch,*cst,*shp,t0,tf);}
   if(s=="CD"      ){ti=new explicitUSL<!remap>                        (*pch,*cst,*shp,t0,tf);}
   if(s=="RKmom"   ){ti=new RK4        < remap>                        (*pch,*cst,*shp,t0,tf);}
   if(s=="RK"      ){ti=new RK4        <!remap>                        (*pch,*cst,*shp,t0,tf);}
   if(s=="linCG"   ){ti=new implicit<conjGrad<impRes> >                (*pch,*cst,*shp,t0,tf,!stat, sym,!remap);}
   if(s=="linGM"   ){ti=new implicit<GMRES<impRes> >                   (*pch,*cst,*shp,t0,tf,!stat,!sym,!remap);}
   if(s=="quasiCG" ){ti=new implicit<NewtonKrylov<conjGrad<NewtSys> > >(*pch,*cst,*shp,t0,tf, stat, sym,!remap);}
   if(s=="NewtCG"  ){ti=new implicit<NewtonKrylov<conjGrad<NewtSys> > >(*pch,*cst,*shp,t0,tf,!stat, sym,!remap);}
   if(s=="asymCG"  ){ti=new implicit<NewtonKrylov<conjGrad<NewtSys> > >(*pch,*cst,*shp,t0,tf,!stat,!sym,!remap);}
   if(s=="quasi"   ){ti=new implicit<NewtonKrylov<GMRES<NewtSys> > >   (*pch,*cst,*shp,t0,tf, stat, sym,!remap);}
   if(s=="dynGM"   ){ti=new implicit<NewtonKrylov<GMRES<NewtSys> > >   (*pch,*cst,*shp,t0,tf,!stat,!sym,!remap);}
   if(s=="momGM"   ){ti=new implicit<NewtonKrylov<GMRES<NewtSys> > >   (*pch,*cst,*shp,t0,tf,!stat, sym, remap);}
   if(s=="dynBCG"  ){ti=new implicit<NewtonKrylov<BiCGstab<NewtSys> > >(*pch,*cst,*shp,t0,tf,!stat, sym,!remap);}
   if(s=="quasiBCG"){ti=new implicit<NewtonKrylov<BiCGstab<NewtSys> > >(*pch,*cst,*shp,t0,tf, stat, sym,!remap);}
   if(ti==NULL)throw stop("makeTimeInt: no time integrator assigned");
}



















