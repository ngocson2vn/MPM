// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include <fstream>
#include "pullinEB.h"
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

struct myPatch:public patch{
   mutable partArray<vec3>surfaceP,posIncP;
   geomStepSeries geoFac;
   timeIntPtr&tmi;
   const shapeSC&shp;
   const realT beamL,beamH,beamB,Ymod,dens,pois,pressure,gravity,permittivity,finalVolt,gap0;
   pullinEB peb;
   mutable realT currVolt;
   ofstream hf,tf;
   int frameCount;
   realT defGivenPres(const realT fperA){                    // fixed-fixed beam with uniform load
      const realT momI=beamB*beamH*beamH*beamH/12.;
      const realT w=-fperA*beamB;                            // load per length
      const realT x=.5*beamL;                                // midpoint of beam
      return w*x*x*(beamL-x)*(beamL-x)/(24.*Ymod*momI);
   }
   realT presGivenFracDef(const realT fracDef){              // fixed-fixed beam with uniform load
      const realT momI=beamB*beamH*beamH*beamH/12.;
      const realT deflect=fracDef*beamH/beamB;               // find load required to deflect by this amount
      const realT x=.5*beamL;                                // midpoint of beam
      return 24.*Ymod*momI*deflect/(x*x*(beamL-x)*(beamL-x));
   }
   realT presGivenGap(const realT gap_mm,const realT volt)const{  // fixed-fixed beam with electro-static load
      const realT gap=1e-3*gap_mm;                           // (meter)
      const realT fperA=permittivity*volt*volt/(2.*gap*gap); // (Newton/meter^2)
      return 1e-6*fperA;                                     // (Newton/mm)
   }
   realT thetaFrac(const realT delt)const{
      const realT oldFrac=geoFac.factorI(tmi->incCount()-2);
      const realT newFrac=geoFac.factorI(tmi->incCount()-1);
      const realT chgFrac=newFrac-oldFrac;
      return oldFrac+(delt/tmi->timeStep())*chgFrac;
   }
   void history();
   bool afterStep();
   void gridVelocityBC(nodeArray<vec3>&,const realT&)const;
   void makeExternalForces(nodeArray<vec3>&,const realT&)const;
   myPatch(timeIntPtr&,constitPtr&,shapePtr&,const realT,const realT,const realT);
   ~myPatch();
};

// This fulfills the virtual function requirement for the caseDesignSC super class.
void myPatch::gridVelocityBC(nodeArray<vec3>&v,const realT&)const{
   for(indexT i=0;i<Nnode();++i){
      //v[i].z=0.; // fouls up lateral stresses - don't do it!
      if(shp.gridPosN(i).x<machEps || shp.gridPosN(i).x>oneEps*beamL){
         v[i]=0.;
      }
   }
}

void myPatch::makeExternalForces(nodeArray<vec3>&gravityN,const realT&delt)const{ // electro pullin
   currVolt=finalVolt*thetaFrac(delt);
   report.param("voltage attempt",currVolt);
   shp.interpolate(posIncP,velN);
   realT mingap=huge;
   const vec3 snapNorm(0,-1,0);
   for(indexT i=0;i<Npart();++i)if(surfNormP[i].y<-machEps){
      const realT refArea=surfAreaP[i].inner(snapNorm);
      const vec3 trialPosP=curPosP[i]+posIncP[i]*tmi->timeStep();
      const realT gap=gap0-abs(trialPosP.y-refPosP[i].y);
      if(gap<mingap)mingap=gap;
      if(gap<gap0*machEps)throw stop("Pull-in detected; exiting early");
      const realT mypressure=presGivenGap(gap,currVolt);
      surfaceP[i]=(mypressure*det(defGradP[i])*refArea*snapNorm.inner(inv(defGradP[i])));
   }
   shp.integrate(surfaceP,gravityN);
   for(indexT i=0;i<Nnode();++i)gravityN[i]/=massN[i];
   //report.param("current gap",mingap);
}

void myPatch::history(){
   if(tmi->incCount()%1==0){
      frameCount+=1;
      hf<<'#'<<frameCount<<'\n';
      //for(indexT i=0;i<Npart();++i){hf<<refPosP[i].x<<'\t'<<velP[i].y<<'\t'<<load*ex->modeShape(mode,refPosP[i].x)<<'\n';}
      for(indexT i=0;i<Npart();++i){hf<<refPosP[i]<<tab<<stressP[i]<<nwl;}
      //for(indexT i=0;i<Npart();++i){hf<<refPosP[i].x<<'\t'<<volumePS[i].xx<<'\t'<<volumePS[i].yy<<'\n';}
      //for(indexT i=0;i<Nnode();++i){hf<<curPosN[i]<<tab<<fintN[i]<<nwl;}
      //for(indexT i=0;i<Nnode();++i)if(abs(curPosN[i].y)<machEps){hf<<curPosN[i].x<<tab<<fintN[i].y<<tab<<fintN[i].x<<nwl;}
      hf<<'\n'<<endl;
   }
}

bool myPatch::afterStep(){
   report.param("voltage success",currVolt);
   realT disp=0.;
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1)disp=curPosP[i].y-refPosP[i].y;
   peb.profileVolt(currVolt);
   tf<<currVolt<<tab<<gap0-abs(disp)<<tab<<gap0-abs(peb.defPos(.5*beamL))<<nwl;
   history();
   return tmi->elapsedTime() < tmi->finalTime;
   //return false;
}

myPatch::~myPatch(){
   report.param("frameCount",frameCount);
   for(indexT i=0;i<Npart();++i)if(probeP[i]==1){
      report.param("finalGap",gap0-abs(curPosP[i].y-refPosP[i].y));
   }
   hf.close(); // various bugs if files not manually closed here
   tf.close();
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
                    pressure(presGivenFracDef(comLineArg("fracDef",.05))),
                    gravity(pressure/(dens*beamH)),
                    permittivity(comLineArg("pmtv",8.854187817e-12) /* (Farad/meter) */ ),
                    finalVolt(comLineArg("volt",100.)),
                    gap0(comLineArg("gap",.00375)),
                    peb(beamL,beamH,beamB,pois,Ymod,gap0,permittivity,101){
   frameCount=0;
   cst=new updatedElastic(Ymod,pois,dens,velGradP,defGradP,volumeP,stressP);
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("quasiCG")),0.,1.);
   report.param("finalTime",tmi->finalTime);
   const vec3 lo(0.,0.,0.);
   const vec3 hi(beamL,beamH,beamB);
   fillIndic(*this,*svar,lo,hi,boxIndic(lo,hi),vec3(2,2,2),false);
   for(indexT i=0;i<Npart();i+=1){
      massP[i]=refVolP[i]*dens;
   }
   partProbe(*svar,*this,vec3(.5*beamL,0.,0.));
   hf.open((string("history.xls")).c_str());  hf.precision(17);
   tf.open((string("tip.xls")).c_str());  tf.precision(17);
   history();
   tf<<0<<tab<<gap0<<tab<<gap0<<nwl;

   ofstream pebfile("pullinEB.xls");  pebfile.precision(14);
   for(realT vo=0.;vo<=finalVolt;vo+=.01*finalVolt){
      peb.profileVolt(vo);
      pebfile<<vo<<tab<<gap0-abs(peb.defPos(.5*beamL))<<endl;
   }
   pebfile.close();
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   realT beamL=comLineArg("beamL",.4);
   realT beamH=comLineArg("beamH",.004);
   realT beamB=comLineArg("beamB",.12);
   int Nthick=comLineArg("Nthick",1);
   const realT cellThick=beamH/realT(Nthick);
   int Nlong=int(round(beamL/cellThick));
   const vec3 lo(0.,-2.*beamH,0.);
   const vec3 hi(beamL,3.*beamH,beamB);
   makeShape(shp,Nlong,5*Nthick,1,lo,hi,comLineArg("shape",string("MPM")));
   myPatch*my=new myPatch(tmi,cst,shp,beamL,beamH,beamB);
   pch=my;
   my->geoFac.init(comLineArg("lastFac",.1),tmi->intEstSteps());
}












