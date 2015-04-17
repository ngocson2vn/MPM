// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

// This problem is taken from Wallstedt and Guilkey
// "An evaluation of explicit time integration schemes for use with the
// generalized interpolation material point method" Journal of
// Computational Physics 2008.

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

struct exactAlign{
   const realT Lx,Ly,Lz,a,b,c,H,rho,lambda,mu;
   exactAlign(const realT&lengthX,
              const realT&lengthY,
              const realT&lengthZ,
              const realT&periodX,
              const realT&periodY,
              const realT&periodZ,
              const realT&amplitude,
              const realT&Ymod,
              const realT&Pois,
              const realT&Dens):
                 Lx(lengthX),
                 Ly(lengthY),
                 Lz(lengthZ),
                 a(periodX),
                 b(periodY),
                 c(periodZ),
                 H(amplitude),
                 rho(Dens),
                 lambda(Ymod*Pois/((1.+Pois)*(1.-2.*Pois))),
                 mu(.5*Ymod/(1.+Pois)){}
   vec3 disp(const vec3&pos,const realT&t){
      const realT&X=pos.x;
      const realT&Y=pos.y;
      const realT&Z=pos.z;
      const realT t2 = 1. / Lx;
      const realT t4 = sin(pi * X * t2);
      const realT t9 = sin(a * pi * t * t2);
      const realT t12 = 1. / Ly;
      const realT t14 = sin(pi * Y * t12);
      const realT t19 = sin(b * pi * t * t12);
      const realT t22 = 1. / Lz;
      const realT t24 = sin(pi * Z * t22);
      const realT t29 = sin(c * pi * t * t22);
      return vec3(H * t4 * t9,
                  H * t14 * t19,
                  H * t24 * t29);
   }
   vec3 vel(const vec3&pos,const realT&t){
      const realT&X=pos.x;
      const realT&Y=pos.y;
      const realT&Z=pos.z;
      const realT t2 = 1. / Lx;
      const realT t4 = sin(pi * X * t2);
      const realT t6 = a * pi;
      const realT t9 = cos(t6 * t * t2);
      const realT t14 = 1. / Ly;
      const realT t16 = sin(pi * Y * t14);
      const realT t18 = b * pi;
      const realT t21 = cos(t18 * t * t14);
      const realT t26 = 1. / Lz;
      const realT t28 = sin(pi * Z * t26);
      const realT t30 = c * pi;
      const realT t33 = cos(t30 * t * t26);
      return vec3(H * t4 * t9 * t6 * t2,
                  H * t16 * t21 * t18 * t14,
                  H * t28 * t33 * t30 * t26);
   }
   vec3 accel(const vec3&pos,const realT&t){
      const realT&X=pos.x;
      const realT&Y=pos.y;
      const realT&Z=pos.z;
      const realT t2 = 1. / Lx;
      const realT t4 = sin(pi * X * t2);
      const realT t9 = sin(a * pi * t * t2);
      const realT t11 = a * a;
      const realT t12 = pi * pi;
      const realT t14 = Lx * Lx;
      const realT t19 = 1. / Ly;
      const realT t21 = sin(pi * Y * t19);
      const realT t26 = sin(b * pi * t * t19);
      const realT t28 = b * b;
      const realT t30 = Ly * Ly;
      const realT t35 = 1. / Lz;
      const realT t37 = sin(pi * Z * t35);
      const realT t42 = sin(c * pi * t * t35);
      const realT t44 = c * c;
      const realT t46 = Lz * Lz;
      return vec3(-H * t4 * t9 * t11 * t12 / t14,
                  -H * t21 * t26 * t28 * t12 / t30,
                  -H * t37 * t42 * t44 * t12 / t46);
   }
   vec3 body(const vec3&pos,const realT&t){
      const realT&X=pos.x;
      const realT&Y=pos.y;
      const realT&Z=pos.z;
      const realT t2 = 1. / Lx;
      const realT t3 = pi * X * t2;
      const realT t4 = sin(t3);
      const realT t5 = H * t4;
      const realT t9 = sin(a * pi * t * t2);
      const realT t11 = a * a;
      const realT t12 = pi * pi;
      const realT t14 = Lx * Lx;
      const realT t18 = 1. / rho;
      const realT t19 = lambda * H;
      const realT t23 = cos(t3);
      const realT t24 = H * t23;
      const realT t28 = 1. + t24 * pi * t2 * t9;
      const realT t32 = Lx + t24 * pi * t9;
      const realT t33 = 1. / t32;
      const realT t38 = 1. / Ly;
      const realT t39 = pi * Y * t38;
      const realT t40 = cos(t39);
      const realT t41 = H * t40;
      const realT t46 = sin(b * pi * t * t38);
      const realT t49 = 1. + t41 * pi * t38 * t46;
      const realT t52 = 1. / Lz;
      const realT t53 = pi * Z * t52;
      const realT t54 = cos(t53);
      const realT t55 = H * t54;
      const realT t60 = sin(c * pi * t * t52);
      const realT t63 = 1. + t55 * pi * t52 * t60;
      const realT t65 = log(t28 * t49 * t63);
      const realT t66 = lambda * t65;
      const realT t67 = t32 * t32;
      const realT t68 = 1. / t67;
      const realT t71 = t5 * t12 * t9;
      const realT t74 = t28 * t28;
      const realT t86 = sin(t39);
      const realT t87 = H * t86;
      const realT t89 = b * b;
      const realT t91 = Ly * Ly;
      const realT t101 = Ly + t41 * pi * t46;
      const realT t102 = 1. / t101;
      const realT t106 = t101 * t101;
      const realT t107 = 1. / t106;
      const realT t110 = t87 * t12 * t46;
      const realT t113 = t49 * t49;
      const realT t125 = sin(t53);
      const realT t126 = H * t125;
      const realT t128 = c * c;
      const realT t130 = Lz * Lz;
      const realT t140 = Lz + t55 * pi * t60;
      const realT t141 = 1. / t140;
      const realT t145 = t140 * t140;
      const realT t146 = 1. / t145;
      const realT t149 = t126 * t12 * t60;
      const realT t152 = t63 * t63;
      return vec3(-t5 * t9 * t11 * t12 / t14 - t18 * (-t19 * t4 * t12 * t2 * t9 / t28 * t33 + t66 * t68 * t71 + mu * t68 * (t74 - 1.) * t71 - 2. * mu * t2 * t33 * t28 * t71),
                  -t87 * t46 * t89 * t12 / t91 - t18 * (-t19 * t86 * t12 * t38 * t46 / t49 * t102 + t66 * t107 * t110 + mu * t107 * (t113 - 1.) * t110 - 2. * mu * t38 * t102 * t49 * t110),
                  -t126 * t60 * t128 * t12 / t130 - t18 * (-t19 * t125 * t12 * t52 * t60 / t63 * t141 + t66 * t146 * t149 + mu * t146 * (t152 - 1.) * t149 - 2. * mu * t52 * t141 * t63 * t149));
   }
   mat3 stress(const vec3&pos,const realT&t){
      const realT&X=pos.x;
      const realT&Y=pos.y;
      const realT&Z=pos.z;
      const realT t2 = 1. / Lx;
      const realT t4 = cos(pi * X * t2);
      const realT t5 = H * t4;
      const realT t10 = sin(a * pi * t * t2);
      const realT t13 = 1. + t5 * pi * t2 * t10;
      const realT t15 = 1. / Ly;
      const realT t17 = cos(pi * Y * t15);
      const realT t18 = H * t17;
      const realT t23 = sin(b * pi * t * t15);
      const realT t26 = 1. + t18 * pi * t15 * t23;
      const realT t29 = 1. / Lz;
      const realT t31 = cos(pi * Z * t29);
      const realT t32 = H * t31;
      const realT t37 = sin(c * pi * t * t29);
      const realT t40 = 1. + t32 * pi * t29 * t37;
      const realT t42 = log(t13 * t26 * t40);
      const realT t43 = lambda * t42;
      const realT t47 = 1. / (Lx + t5 * pi * t10);
      const realT t51 = t13 * t13;
      const realT t59 = 1. / (Ly + t18 * pi * t23);
      const realT t63 = t26 * t26;
      const realT t71 = 1. / (Lz + t32 * pi * t37);
      const realT t75 = t40 * t40;
      return mat3(t43 * Lx * t47 + mu * Lx * t47 * (t51 - 1.),0.,0.,
                  0.,t43 * Ly * t59 + mu * Ly * t59 * (t63 - 1.),0.,
                  0.,0.,t43 * Lz * t71 + mu * Lz * t71 * (t75 - 1.));
   }
};

// case definition
// all problem-specific information is pulled under this umbrella
struct myPatch:public patch{
   mutable partArray<vec3>sourceP;
   ofstream hf;
   realT Li,L2,L1,Ymod,dens,pois,load,vwav;
   int frameCount;
   neoHookean*nh;
   exactAlign*ex;
   timeIntPtr&tmi;
   shapePtr&shp;
   bool single;
   bool afterStep();
   void gridVelocityBC(nodeArray<vec3>&,const realT&)const;
   void makeExternalForces(nodeArray<vec3>&,const realT&)const;
   myPatch(timeIntPtr&tmi,constitPtr&cst,shapePtr&s);
   ~myPatch();
};

// This fulfills the virtual function requirement for the caseDesignSC super class.
void myPatch::makeExternalForces(nodeArray<vec3>&gravityN,const realT&)const{
   for(indexT i=0;i<Npart();++i)sourceP[i]=massP[i]*ex->body(refPosP[i],tmi->elapsedTime()-tmi->timeStep());
   shp->integrate(sourceP,gravityN);
   for(indexT i=0;i<Nnode();++i)gravityN[i]/=massN[i];
}

// This fulfills the virtual function requirement for the caseDesignSC super class.
void myPatch::gridVelocityBC(nodeArray<vec3>&v,const realT&)const{
   for(indexT i=0;i<Nnode();++i){
      if(shp->gridPosN(i).x<machEps||shp->gridPosN(i).x>1.-machEps)v[i].x=0.;
      if(shp->gridPosN(i).y<machEps||shp->gridPosN(i).y>1.-machEps)v[i].y=0.;
      if(shp->gridPosN(i).z<machEps||shp->gridPosN(i).z>1.-machEps)v[i].z=0.;
   }
}

// This fulfills the virtual function requirement for the caseDesignSC super class.
bool myPatch::afterStep(){
   // If its a single run, we write data to disk.  Otherwise we
   // only collect measures of error.
   if(single){
      frameCount+=1;
      hf<<'#'<<frameCount<<'\n';
      for(indexT i=0;i<Npart();++i){hf<<curPosP[i]<<'\t'<<det(defGradP[i])<<'\n';}
      //for(int i=0;i<Nnode();++i){hf<<pch->curPosN[i]<<'\t'<<mag(pch->ga[i])<<'\n';}
      hf<<"\n\n";
   }
   // Here we measure error as the magnitude of the vector of
   // displacement error: ||u_actual - u_exact||
   realT Lisum=0.;
   realT L2sum=0.;
   realT L1sum=0.;
   for(indexT i=0;i<Npart();++i){
      const realT err=radius(curPosP[i]-refPosP[i],ex->disp(refPosP[i],tmi->elapsedTime()));
      L1sum+=abs(err);
      L2sum+=err*err;
      Lisum=(err>Lisum?err:Lisum);
   }
   realT totP=realT(Npart());
   MPI_Allreduce(MPI_IN_PLACE,&totP ,1,MPI_DOUBLE,MPI_SUM,cartComm);
   MPI_Allreduce(MPI_IN_PLACE,&L1sum,1,MPI_DOUBLE,MPI_SUM,cartComm);
   MPI_Allreduce(MPI_IN_PLACE,&L2sum,1,MPI_DOUBLE,MPI_SUM,cartComm);
   MPI_Allreduce(MPI_IN_PLACE,&Lisum,1,MPI_DOUBLE,MPI_MAX,cartComm);
   L1sum=L1sum/totP;
   L2sum=sqrt(L2sum/totP);
   // continually look for the highest error
   L1=(L1sum>L1?L1sum:L1);  // L-1 norm
   L2=(L2sum>L2?L2sum:L2);  // L-2 norm
   Li=(Lisum>Li?Lisum:Li);  // L-infinity norm
   return tmi->elapsedTime() < tmi->finalTime;
   //return false;
}

// This fulfills the virtual function requirement for the caseDesignSC super class.
myPatch::~myPatch(){
   report.param("frameCount",frameCount);
   report.param("L1norm",L1);
   report.param("L2norm",L2);
   report.param("Linfnorm",Li);
   if(cartRank==0)cerr<<"L1-norm   "<<L1<<endl;
   if(cartRank==0)cerr<<"Linf-norm "<<Li<<endl;
   hf.close(); // various bugs if files not manually closed here
}

// This is the myPatch constructor
myPatch::myPatch(timeIntPtr&tm,constitPtr&cst,shapePtr&svar):tmi(tm),shp(svar){
   const vec3 lo(0.,0.,0.);
   const vec3 hi(1.,1.,1.);
   int Nc=comLineArg("Ncell",20);
   makeShape(shp,Nc,Nc,Nc,lo,hi,comLineArg("shape",string("GIMP")));

   Li=L2=L1=0.;
   single=true;
   frameCount=0;

   load=comLineArg("load",.1);
   Ymod=1e4;
   dens=1.;
   vwav=sqrt(Ymod/dens);
   pois=.3;

   // make a new exact solution class, if one is defined
   //ex=new exactAlign(1.,1.,1.,vwav,vwav,vwav,load,Ymod,pois,dens);
   ex=new exactAlign(1.,1.,1.,vwav,vwav,vwav,load,Ymod,pois,dens);

   // This function fills a rectangle with particles.
   boxIndic indic(shp->regionLo,shp->regionHi);
   fillIndic(*this,*svar,shp->regionLo,shp->regionHi,indic,vec3(2,2,2),false);

   // initialize constitutive model
   cst=nh=new neoHookean(Ymod,pois,dens,velGradP,defGradP,volumeP,refVolP,stressP);

   // We initialize the simulation time to zero.
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("CD")),0.,2./vwav);

   // Initializing particles to the exact values of the manufactured solution
   // if there is one, or to some known initial state
   for(indexT i=0;i<Npart();i+=1){
      massP[i]=refVolP[i]*dens;
      velP[i]=ex->vel(refPosP[i],tmi->elapsedTime());
      //velP[i]=ex->vel(refPosP[i],tmi->elapsedTime()-.5*tmi->timeStep()); // not necessary?
   }
   report.param("Ymod",Ymod);
   report.param("dens",dens);
   report.param("pois",pois);
   report.param("load",load);
   report.param("vwav",vwav);
   hf.open("history.xls");  hf.precision(14);
   single=comLineArg("single",true);
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){pch=new myPatch(tmi,cst,shp);}


