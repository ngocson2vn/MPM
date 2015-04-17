// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef TIMEINT_H
#define TIMEINT_H

#include<fstream>
#include"main.h"
using namespace std;

inline timeIntSC::timeIntSC(const patch&,const constitutiveSC&cst,const shapeSC&shp,realT&t0,realT&tf):initialTime(t0),finalTime(tf){
   _incCount=0;
   _elapsedTime=initialTime;
   const realT CFL=comLineArg("CFL",.5);
   realT minh=shp.cell.x;
   if(shp.cell.y<minh)minh=shp.cell.y;
   if(shp.cell.z<minh)minh=shp.cell.z;
   nominalStep=CFL*minh/cst.waveSpeed();
   const realT totTime=finalTime-initialTime;
   Nstep=comLineArg("Nstep",totTime/nominalStep);     // allow Nstep to modify nominalStep
   nominalStep=comLineArg("dtOveride",totTime/Nstep); // allow dtOveride to modify nominalStep
   Nstep=totTime/nominalStep;                         // in case nominal step was modified
   // If nominalStep is not modified, then nominalStep = totTime/(totTime/nominalStep) = nominalStep
   report.param("timeStep",nominalStep);
}

inline void timeIntSC::nextStep(const patch&){
   const realT allowedRoundoff=nominalStep*machEps*realT(_incCount);
   const realT timeRemaining=finalTime-_elapsedTime;
   _timeStep=(timeRemaining<nominalStep+allowedRoundoff?timeRemaining:nominalStep);
   _elapsedTime+=_timeStep;
   ++_incCount;
}

class timeIntOps{
   partArray<vec3>incP,momP;
   nodeArray<vec3>momN;
   realT minMass;
   shapeSC&shp;
public:
   timeIntOps(shapeSC&s):shp(s){}
   void setMinMass(const patch&pch);
   void integrateMass(const partArray<realT>&massP,nodeArray<realT>&massN);
   void mapVelPtoVelN(const patch&pch,const realT&dt,const partArray<vec3>&velP,nodeArray<vec3>&velN);
   void solveMomentum(const patch&pch,
                      const realT&dt,
                      nodeArray<vec3>&fintN,
                      nodeArray<vec3>&velN,
                      nodeArray<vec3>&velIncN);
   void advancePosition(const nodeArray<vec3>&velN,const realT&dt,partArray<vec3>&curPosP);
   void advanceVelocity(const nodeArray<vec3>&velIncN,const realT&dt,partArray<vec3>&velP);
   void advanceConstit(patch&pch,constitutiveSC&cst,const realT&dt);
};

//////// implicit methods /////////////////////

class impRes:public timeIntOps{
   partArray<vec3>velCopyP;
   nodeArray<vec3>vTmp,rhs;
   realT bnrm2;
public:
   const realT theta;
private:
   patch&pch; // const?
   constitutiveSC&cst;
   shapeSC&shp; // const?
   const nodeArray<vec3>&fintStartN;
   const nodeArray<vec3>&vStart;
   indexT countRes;
   const bool stat,sym,remap;
   void common(nodeArray<vec3>&v,nodeArray<vec3>&res,const realT&delt){
      pch.gridVelocityBC(v,delt);
      for(indexT i=0;i<Nnode();++i)vTmp[i]=v[i];
      if(remap){
         for(indexT i=0;i<Nnode();++i)v[i]=(v[i]-vStart[i]); // theta?
         for(indexT i=0;i<Npart();++i)velCopyP[i]=pch.velP[i];
         advanceVelocity(v,0.,velCopyP);
         mapVelPtoVelN(pch,0.,velCopyP,v);
      }
      cst.revert();
      shp.gradient(pch.velGradP,v);
      cst.update(delt);
      shp.divergence(pch.stressP,pch.volumeP,pch.fintN);
      for(indexT i=0;i<Nnode();++i)res[i]=-delt*(pch.fintN[i]-fintStartN[i])/pch.massN[i];
      if(!stat){
         if(remap)for(indexT i=0;i<Nnode();++i)res[i]+=vTmp[i];
         else for(indexT i=0;i<Nnode();++i)res[i]+=v[i];
      }
      pch.gridVelocityBC(res,delt);
      ++countRes;
   }
public:
   impRes(patch&p,
          constitutiveSC&c,
          shapeSC&s,
          const nodeArray<vec3>&v0,
          const nodeArray<vec3>&f0,
          const bool st,
          const bool sy,
          const bool re):timeIntOps(s),
                         theta(comLineArg("theta",(st?1.:.5))),
                         pch(p),
                         cst(c),
                         shp(s),
                         fintStartN(f0),
                         vStart(v0),
                         stat(st),
                         sym(sy),
                         remap(re){bnrm2=0.;}
   realT normb(){return bnrm2;}
   void makeLHS(nodeArray<vec3>&v,nodeArray<vec3>&res,const realT&delt){
      common(v,res,delt);
      if(sym)for(indexT i=0;i<Nnode();++i)res[i]*=pch.massN[i]/bnrm2;
   }
   void makeRHSb(nodeArray<vec3>&res,const realT&delt){
      for(indexT i=0;i<Nnode();++i)res[i]=-delt*(pch.bodyAccelN[i] + fintStartN[i]/pch.massN[i]);
      if(!stat)for(indexT i=0;i<Nnode();++i)res[i]-=vStart[i];
      for(indexT i=0;i<Nnode();++i)res[i]=-res[i];
      pch.gridVelocityBC(res,delt);
   }
   void makeRes(nodeArray<vec3>&v,nodeArray<vec3>&res,const realT&delt){
      common(v,res,delt);
      makeRHSb(rhs,delt);
      //bnrm2=shp.dot(rhs,rhs);
      bnrm2=1.;
      for(indexT i=0;i<Nnode();++i)res[i]=rhs[i]-res[i];
      if(sym)for(indexT i=0;i<Nnode();++i)res[i]*=pch.massN[i]/bnrm2;
   }
};

class NewtSys:public impRes{
   nodeArray<vec3>vPert,Zpert,Zcurr;
   nodeArray<vec3>&velN;
   const shapeSC&shp;
public:
   NewtSys(patch&p,
           constitutiveSC&c,
           shapeSC&sh,
           const nodeArray<vec3>&v0,
           const nodeArray<vec3>&f0,
           const bool st,
           const bool sy,
           const bool re):
              impRes(p,c,sh,v0,f0,st,sy,re),
              velN(p.velN),
              shp(sh){}
   void NewtonRes(nodeArray<vec3>&v,nodeArray<vec3>&r,const realT&delt){impRes::makeRes(v,r,delt);}
   void makeLHS(const nodeArray<vec3>&s,nodeArray<vec3>&DZ,const realT&delt){
      const realT dotss=shp.dot(s,s);
      const realT dotvv=shp.dot(velN,velN);
      const realT fd=sqrt(machEps) + sqrt((1.+dotvv)*machEps/dotss);
      assert(fd==fd);
      impRes::makeRes(velN,Zcurr,delt);
      for(indexT i=0;i<Nnode();++i)vPert[i]=velN[i]+fd*s[i];
      impRes::makeRes(vPert,Zpert,delt);
      for(indexT i=0;i<Nnode();++i)DZ[i]=(Zpert[i]-Zcurr[i])/fd;
   }
   void makeRHSb(nodeArray<vec3>&Z,const realT&delt){ // first argument ignored
      impRes::makeRes(velN,Z,delt);
   }
   void makeRes(nodeArray<vec3>&,nodeArray<vec3>&Z,const realT&delt){ // first argument ignored
      impRes::makeRes(velN,Z,delt); // not valid for GMRES re-starting where full res=-F-Js is needed
      for(indexT i=0;i<Nnode();++i)Z[i]=-Z[i];
   }
};

//////// Krylov solvers ///////////////////

template<typename sysT>
class BiCGstab:public sysT{
   nodeArray<vec3>r,r_tld,P,V,S,T;
   const realT linTol;
   const indexT linIter;
   const shapeSC&shp;
   indexT Kiter;
public:
   BiCGstab(patch&p,
            constitutiveSC&c,
            shapeSC&sh,
            const nodeArray<vec3>&v0,
            const nodeArray<vec3>&f0,
            const bool st,
            const bool sy,
            const bool re):
               sysT(p,c,sh,v0,f0,st,sy,re),
               linTol(comLineArg("linTol",1e-8)),
               linIter(comLineArg("linIter",256)),
               shp(sh){}
   indexT solve(nodeArray<vec3>&primary,const realT&delt){
      sysT::makeRes(primary,r,delt);
      if(sysT::normb() == 0.)return 0;
      if(shp.dot(r,r) / sysT::normb() < linTol)return 0;
      realT omega = 1.,rho=huge,rho_1=huge,alpha=huge;
      for(indexT i=0;i<Nnode();++i)r_tld[i]=r[i];
      for(Kiter=0;Kiter<linIter;++Kiter){
         if(Kiter%(1+indexT(realT(1000000)/realT(Nnode())))==0){
            report.param("current Kiter",Kiter);
            report.write();
         }
         rho=shp.dot(r_tld,r);
         if(rho == 0.)break;
         if(Kiter != 0){
            realT beta=(rho/rho_1)*(alpha/omega);
            for(indexT i=0;i<Nnode();++i)P[i]=r[i]+beta*(P[i]-omega*V[i]);
         }
         else{
            for(indexT i=0;i<Nnode();++i)P[i]=r[i];
         }
         sysT::makeLHS(P,V,delt);
         alpha=rho/shp.dot(r_tld,V);
         for(indexT i=0;i<Nnode();++i)S[i]=r[i]-alpha*V[i];
         if(shp.dot(S,S) < linTol){
            for(indexT i=0;i<Nnode();++i)primary[i]+=alpha*P[i];
            break;
         }
         sysT::makeLHS(S,T,delt);
         const realT dotTS=shp.dot(T,S);
         const realT dotTT=shp.dot(T,T);
         omega=dotTS/dotTT;
         for(indexT i=0;i<Nnode();++i)primary[i]+=alpha*P[i]+omega*S[i];
         for(indexT i=0;i<Nnode();++i)r[i]=S[i]-omega*T[i];
         if(shp.dot(r,r)/sysT::normb() <= linTol)break;
         if(omega == 0.)break;
         rho_1=rho;
      }
      return Kiter;
   }
};

template<typename sysT>
class conjGrad:public sysT{
   nodeArray<vec3>res,P,Q;
   map<indexT,indexT>histogram;
   const realT linTol;
   const indexT linIter;
   const shapeSC&shp;
   indexT Kiter,totIter;
public:
   conjGrad(patch&p,
            constitutiveSC&c,
            shapeSC&sh,
            const nodeArray<vec3>&v0,
            const nodeArray<vec3>&f0,
            const bool st,
            const bool sy,
            const bool re):
               sysT(p,c,sh,v0,f0,st,sy,re),
               linTol(comLineArg("linTol",1e-8)),
               linIter(comLineArg("linIter",256)),
               shp(sh){totIter=0;}

   ~conjGrad(){
      cerr<<"Histogram of required CG iterations per call"<<endl;
      for(indexT i=0;i<histogram.size();++i){
         if(histogram[i]>0)cerr<<"  iters "<<setw(5)<<i+1<<"  times reached "<<histogram[i]<<endl;
      }
      report.param("total CG iterations",totIter);
      cerr<<"total CG iterations "<<totIter<<endl;
   }

   indexT solve(nodeArray<vec3>&primary,const realT&delt){
      sysT::makeRes(primary,res,delt);
      P=res;
      realT rhoNew2=shp.dot(res,res);
      const realT rhoNewton=sqrt(rhoNew2);
      for(Kiter=0;Kiter<linIter;++Kiter){
         ++totIter;
         if(Kiter%(1+indexT(realT(1000000)/realT(Nnode())))==0){
            report.param("current Kiter",Kiter);
            report.write();
         }
         sysT::makeLHS(P,Q,delt);
         const realT dotpq=shp.dot(P,Q);
         const realT alpha=rhoNew2/dotpq;
         for(indexT i=0;i<Nnode();++i)primary[i]+=alpha*P[i];
         for(indexT i=0;i<Nnode();++i)res[i]-=alpha*Q[i];
         const realT rhoOld2=rhoNew2;
         rhoNew2=shp.dot(res,res);
         const realT ratio=sqrt(rhoNew2)/rhoNewton;
         if(ratio<linTol)break; // achieved the goal
         const realT beta=rhoNew2/rhoOld2;
         for(indexT i=0;i<Nnode();++i)P[i]=res[i]+beta*P[i];
      }
      histogram[Kiter-1]+=1; // hope the initial value is 0
      return Kiter;
   }
};

template<typename sysT>
class GMRES:public sysT{
   vector<realT>gNorm,cGiv,sGiv,yMin;
   vector<vector<realT> >Hess;
   vector<nodeArray<vec3> >qRes;
   vector<indexT>histogram;
   const realT linTol;
   const indexT linIter;
   const shapeSC&shp;
   indexT countBasis,totIter;
   nodeArray<vec3>res;
   void multGivens(const indexT k,realT&x,realT&y){
      const realT tmp=-sGiv[k]*x+cGiv[k]*y;
      x=cGiv[k]*x+sGiv[k]*y;
      y=tmp;
   }
public:
   GMRES(patch&p,
         constitutiveSC&c,
         shapeSC&sh,
         const nodeArray<vec3>&v0,
         const nodeArray<vec3>&f0,
         const bool st,
         const bool sy,
         const bool re):
            sysT(p,c,sh,v0,f0,st,sy,re),
            linTol(comLineArg("linTol",1e-8)),
            linIter(comLineArg("linIter",128)),
            shp(sh){
      report.param("linIter",linIter);
      cGiv.resize(linIter);
      sGiv.resize(linIter);
      gNorm.resize(linIter+1);
      yMin.resize(linIter+1);
      qRes.resize(linIter+1,nodeArray<vec3>());
      Hess.resize(linIter+1,vector<realT>(linIter));
      histogram.resize(linIter);
      for(indexT i=0;i<histogram.size();++i)histogram[i]=0;
      countBasis=0;
      totIter=0;
   }
   ~GMRES(){
      cerr<<"Histogram of required GMRES iterations per call"<<endl;
      for(indexT i=0;i<histogram.size();++i){
         if(histogram[i]>0)cerr<<"  iters "<<setw(5)<<i+1<<"  times reached "<<histogram[i]<<endl;
      }
      report.param("countBasis",countBasis);
      report.param("total GMRES iterations",totIter);
      cerr<<"total GMRES iterations "<<totIter<<endl;
   }
   indexT solve(nodeArray<vec3>&primary,const realT&delt){
      sysT::makeRes(primary,res,delt);
      realT rhoCurr=sqrt(shp.dot(res,res));
      const realT rhoNewton=rhoCurr;
      realT rhoPrev=huge;
      for(indexT i=0;i<Nnode();++i)qRes[0][i]=res[i]/rhoCurr;
      for(indexT i=0;i<Hess.size();++i)Hess[i].assign(Hess[i].size(),0.);
      gNorm.assign(gNorm.size(),0.);
      gNorm[0]=rhoCurr;
      int Kiter=-1;
      while(Kiter<int(linIter)-1 && rhoCurr>rhoNewton*linTol){
         ++totIter;
         rhoPrev=rhoCurr;
         if(Kiter%100==0){
            report.param("current Kiter",Kiter);
            report.write();
         }
         ++Kiter;
         sysT::makeLHS(qRes[Kiter],qRes[Kiter+1],delt);
         for(int j=0;j<=Kiter;++j){
            Hess[j][Kiter]=shp.dot(qRes[Kiter+1],qRes[j]);
            for(indexT i=0;i<Nnode();++i)qRes[Kiter+1][i]-=Hess[j][Kiter]*qRes[j][i];
            ++countBasis;
         }
         Hess[Kiter+1][Kiter]=sqrt(shp.dot(qRes[Kiter+1],qRes[Kiter+1]));
         for(indexT i=0;i<Nnode();++i)qRes[Kiter+1][i]/=Hess[Kiter+1][Kiter];
         for(int i=0;i<Kiter;++i)multGivens(i,Hess[i][Kiter],Hess[i+1][Kiter]);
         const realT h0=Hess[Kiter][Kiter];
         const realT h1=Hess[Kiter+1][Kiter];
         const realT nu=sqrt(h0*h0+h1*h1);
         cGiv[Kiter]=Hess[Kiter][Kiter]/nu;
         sGiv[Kiter]=Hess[Kiter+1][Kiter]/nu;
         Hess[Kiter][Kiter]=cGiv[Kiter]*Hess[Kiter][Kiter]+sGiv[Kiter]*Hess[Kiter+1][Kiter];
         Hess[Kiter+1][Kiter]=0.;
         multGivens(Kiter,gNorm[Kiter],gNorm[Kiter+1]);
         rhoCurr=abs(gNorm[Kiter+1]);
         if(Kiter>0 && rhoCurr>rhoPrev){        // GMRES bottomed out - we need to make a choice
            if(rhoCurr>rhoNewton){
               cerr<<"!! Discarding GMRES result; primary=0"<<endl;
               return -1;                       // discard entire GMRES result - primary=0
            }
            --Kiter;                            // discard this but keep the previous
            break;
         }
      }
      ++histogram[Kiter];
      yMin[Kiter]=gNorm[Kiter]/Hess[Kiter][Kiter];
      for(int i=Kiter-1;i>-1;--i){
         yMin[i]=gNorm[i];
         for(int j=i+1;j<=Kiter;++j)yMin[i]-=Hess[i][j]*yMin[j];
         yMin[i]/=Hess[i][i];
      }
      for(indexT i=0;i<Nnode();++i){
         for(int j=0;j<=Kiter;++j)primary[i]+=qRes[j][i]*yMin[j];
      }
      return Kiter;
   }
};

//////// Newton-Krylov ////////////////////////

template<typename solverT>
class NewtonKrylov:public solverT{
   ofstream newtFile;
   nodeArray<vec3>vInc,Zcurr,vTry;
   const realT ratioNewton,ignoreResidualBelow;
   realT rhoStart;
   const indexT NewtLim,maxLoadIter;
   patch&pch;
   constitutiveSC&cst;
   const nodeArray<vec3>&vStart;
   const shapeSC&shp;
public:
   NewtonKrylov(patch&p,
          constitutiveSC&c,
          shapeSC&sh,
          const nodeArray<vec3>&v0,
          const nodeArray<vec3>&f0,
          const bool st,
          const bool sy,
          const bool re):
             solverT(p,c,sh,v0,f0,st,sy,re),
             ratioNewton(comLineArg("NewtTol",1e-6)),
             ignoreResidualBelow(comLineArg("ignore",tiny)),
             rhoStart(-huge),
             NewtLim(comLineArg("NewtLim",8)),
             maxLoadIter(comLineArg("maxLoadIter",1)),
             pch(p),
             cst(c),
             vStart(v0),
             shp(sh){
      newtFile.open("newton.xls");  newtFile.precision(14);
   }
   ~NewtonKrylov(){
      newtFile.close();
   }
   void solve(nodeArray<vec3>&velN2,const realT&delt){
      realT loadResPrev=huge;
      for(indexT Nload=0;Nload<maxLoadIter;++Nload){
         cst.revert();
         pch.makeExternalForces(pch.bodyAccelN,delt);
         solverT::NewtonRes(velN2,Zcurr,delt);
         realT rhoCurr=sqrt(shp.dot(Zcurr,Zcurr));
         realT rhoPrev=huge;
         if(Nload==0){
            rhoStart=rhoCurr;
            report.param("rhoStart",rhoStart);
         }
         else{
            if(rhoCurr<rhoStart*ratioNewton)break;
            if(rhoCurr/loadResPrev>1.)throw stop("nonlinear load convergence failed");
            report.param("nonlinear load loop",Nload);
            report.param("nonlinear load loop tol start",rhoCurr/rhoStart);
            report.param("nonlinear load loop tol prev",rhoCurr/loadResPrev);
         }
         report.write();
         indexT Nnewton=0;
         loadResPrev=rhoCurr;
         while(rhoCurr>ignoreResidualBelow && rhoCurr>=rhoStart*ratioNewton){
            ++Nnewton;
            report.param("Newton step",Nnewton);
            for(indexT i=0;i<Nnode();++i)vInc[i]=0.;
            const indexT Kiter2=solverT::solve(vInc,delt);
            report.param("ending Kiter",Kiter2);
            for(indexT ArmijoCount=0;ArmijoCount<8;++ArmijoCount){
               for(indexT i=0;i<Nnode();++i)vTry[i]=velN2[i]+vInc[i]/pow(2.,ArmijoCount);
               solverT::NewtonRes(vTry,Zcurr,delt);
               rhoPrev=rhoCurr;
               rhoCurr=sqrt(shp.dot(Zcurr,Zcurr));
               newtFile<<rhoCurr/rhoStart<<tab<<rhoCurr/rhoPrev<<tab<<Kiter2+1<<endl;
               report.param("tol start",rhoCurr/rhoStart);
               report.param("tol prev",rhoCurr/rhoPrev);
               report.write();
               if(rhoCurr<rhoStart*ratioNewton){velN2=vTry;break;} // achieved the goal - go to next time step
               if(rhoCurr<rhoPrev)             {velN2=vTry;break;} // made progress - try another Newton
               if(Nnewton>1 && ArmijoCount>=7){                    // Armijo failed - discard result
                  report.param("(Armijo) residual not satisfied",rhoCurr/rhoStart);
                  cerr<<"(Armijo) residual not satisfied "<<rhoCurr/rhoStart<<endl;
                  return;
               }
               report.param("Invoking Armijo tol prev",rhoCurr/rhoPrev);
            }
            if(Nnewton>=NewtLim){
               report.param("(Newton) residual not satisfied",rhoCurr/rhoStart);
               cerr<<"(Newton) residual not satisfied "<<rhoCurr/rhoStart<<endl;
               return;
            }
         }
      }
   }
};

//////// drivers ///////////////////////

template<bool remap>
class RK4:public timeIntSC,public timeIntOps{
   nodeArray<vec3>inc1,inc2,inc3,inc4,vTmp,vStart;
   partArray<vec3>velCopyP;
   realT prevStep;
   patch&pch;
   constitutiveSC&cst;
   shapeSC&shp;
   void stage(const nodeArray<vec3>&prev,nodeArray<vec3>&inc,const realT&coeff){
      const realT&dt=timeStep();
      for(indexT i=0;i<Nnode();++i)vTmp[i]=pch.velN[i]+coeff*prev[i];
      pch.gridVelocityBC(vTmp,coeff*dt);
      cst.revert();
      if(remap){
         for(indexT i=0;i<Nnode();++i)vTmp[i]=(vTmp[i]-vStart[i]);
         for(indexT i=0;i<Npart();++i)velCopyP[i]=pch.velP[i];
         advanceVelocity(vTmp,0.,velCopyP);
         mapVelPtoVelN(pch,0.,velCopyP,vTmp);
      }
      shp.gradient(pch.velGradP,vTmp);
      cst.update(coeff*dt);
      shp.divergence(pch.stressP,pch.volumeP,pch.fintN);
      pch.makeExternalForces(pch.bodyAccelN,coeff*dt);
      for(indexT i=0;i<Nnode();++i)inc[i]=(pch.bodyAccelN[i]+pch.fintN[i]/pch.massN[i])*dt;
   }
public:
   RK4(patch&p,constitutiveSC&c,shapeSC&s,realT&t0,realT&tf):timeIntSC(p,c,s,t0,tf),timeIntOps(s),pch(p),cst(c),shp(s){prevStep=0.;}
   void advance(const shapeSC&){
      const realT&dt=timeStep();
      if(incCount()==1)setMinMass(pch);
      integrateMass(pch.massP,pch.massN);
      mapVelPtoVelN(pch,dt,pch.velP,pch.velN);
      cst.save();
      for(indexT i=0;i<Nnode();++i)vStart[i]=pch.velN[i];
      stage(pch.velN,inc1,0.);
      stage(inc1    ,inc2,.5);
      stage(inc2    ,inc3,.5);
      stage(inc3    ,inc4,1.);
      cst.revert();
      for(indexT i=0;i<Nnode();++i)vTmp[i]=pch.velN[i];
      for(indexT i=0;i<Nnode();++i)pch.velN[i]+=((inc1[i]+2.*inc2[i])+(2.*inc3[i]+inc4[i]))/6.;
      pch.gridVelocityBC(pch.velN,dt);
      for(indexT i=0;i<Nnode();++i)pch.velIncN[i]=(pch.velN[i]-vTmp[i]);
      advancePosition(pch.velN,dt,pch.curPosP);
      advanceVelocity(pch.velIncN,dt,pch.velP);
      if(remap)mapVelPtoVelN(pch,dt,pch.velP,pch.velN);
      advanceConstit(pch,cst,dt);
      prevStep=dt;
   }
};

template<bool remap>
class explicitUSL:public timeIntSC,public timeIntOps{
   realT prevStep;
   patch&pch;
   constitutiveSC&cst;
public:
   explicitUSL(patch&p,constitutiveSC&c,shapeSC&s,realT&t0,realT&tf):timeIntSC(p,c,s,t0,tf),timeIntOps(s),pch(p),cst(c){prevStep=0.;}
   void advance(const shapeSC&){
      const realT&dt=timeStep();
      if(incCount()==1)setMinMass(pch);
      integrateMass(pch.massP,pch.massN);
      mapVelPtoVelN(pch,dt,pch.velP,pch.velN);
      pch.makeExternalForces(pch.bodyAccelN,dt);
      solveMomentum(pch,.5*(prevStep+dt),pch.fintN,pch.velN,pch.velIncN);
      advancePosition(pch.velN,dt,pch.curPosP);
      advanceVelocity(pch.velIncN,dt,pch.velP);
      if(remap)mapVelPtoVelN(pch,dt,pch.velP,pch.velN);
      advanceConstit(pch,cst,dt);
      prevStep=dt;
   }
};

template<typename NewtT>
class implicit:public NewtT,public timeIntSC{
   nodeArray<vec3>vStart,fintStartN;
   patch&pch;
   constitutiveSC&cst;
   const bool stat,remap;
public:
   implicit(patch&p,
                   constitutiveSC&cc,
                   shapeSC&sh,
                   realT&t0,realT&tf,
                   const bool st,
                   const bool sy,
                   const bool re):
                      NewtT(p,cc,sh,vStart,fintStartN,st,sy,re),
                      timeIntSC(p,cc,sh,t0,tf),
                      pch(p),
                      cst(cc),
                      stat(st),
                      remap(re){}
   void advance(const shapeSC&shp2){
      const realT&dt=timeStep();
      if(incCount()==1)NewtT::setMinMass(pch);
      NewtT::integrateMass(pch.massP,pch.massN);
      pch.makeExternalForces(pch.bodyAccelN,NewtT::theta*dt);
      NewtT::mapVelPtoVelN(pch,dt,pch.velP,pch.velN);
      shp2.divergence(pch.stressP,pch.volumeP,fintStartN);
      vStart=pch.velN;
      cst.save();
      NewtT::solve(pch.velN,NewtT::theta*dt);
      cst.revert();
      if(!stat){
         if(remap)for(indexT i=0;i<Nnode();++i)pch.velIncN[i]=(pch.velN[i]-vStart[i]);
         else     for(indexT i=0;i<Nnode();++i)pch.velIncN[i]=(pch.velN[i]-vStart[i])/NewtT::theta;
         NewtT::advanceVelocity(pch.velIncN,dt,pch.velP);
      }
      NewtT::advancePosition(pch.velN,dt,pch.curPosP);
      if(remap)NewtT::mapVelPtoVelN(pch,dt,pch.velP,pch.velN);
      NewtT::advanceConstit(pch,cst,dt);
   }
};


















#endif
