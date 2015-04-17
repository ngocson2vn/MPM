// 2011 Philip Wallstedt

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

struct myPatch:public patch{
   ofstream hf;
   timeIntPtr&tmi;
   const shapeSC&shp;
   bool afterStep();
   myPatch(timeIntPtr&,constitPtr&cst,shapePtr&s);
   //~myPatch(){
      //for(indexT i=0;i<Npart();++i){
         //report.param((toStr("part")+toStr(i)).c_str(),curPosP[i]-refPosP[i]);
      //}
   //}
};

bool myPatch::afterStep(){
   if(cartRank==0)for(indexT i=0;i<Nnode();++i){
      if(massN[i]>.0001)hf<<shp.gridPosN(i)<<tab<<mag(velN[i]*massN[i])<<nwl;
      //cerr<<i<<tab<<massN[i]<<nwl;
   }
   if(cartRank==0)cerr<<nwl<<"step "<<tmi->incCount()<<" ------------"<<endl;
   for(int c=0;c<cartSize;++c){
      MPI_Barrier(cartComm);
      if(c==cartRank)for(indexT i=0;i<Npart();++i)cerr<<cartRank<<tab<<curPosP[i]-refPosP[i]<<nwl;
      cerr<<flush;
      MPI_Barrier(cartComm);
   }
   //cerr<<"part position0 "<<curPosP[0]<<tab<<tmi->elapsedTime()<<endl;
   //return false;
   return tmi->incCount() < 5;
}

myPatch::myPatch(timeIntPtr&tm,constitPtr&cst,shapePtr&svar):tmi(tm),shp(*svar){
   fillIndic(*this,*svar,vec3(),vec3(1,1,1),boxIndic(vec3(),vec3(1,1,1)),vec3(1,1,1),false);
   cst=new neoHookean(1000.,.25,1.,velGradP,defGradP,volumeP,refVolP,stressP);
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("CD")),0.,1./32.);
   const realT load=comLineArg("load",32.);
   for(indexT i=0;i<Npart();++i){
      massP[i]=1.;
      velP[i]=load*vec3(1,1,1);
   }
   //cerr<<"part position0 "<<curPosP[0]<<tab<<tmi->elapsedTime()<<endl;
   if(cartRank==0)hf.open("history.xls");  hf.precision(14);
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   const vec3 lo(0.,0.,0.);
   const vec3 hi(2.,2.,2.);
   int Nc=comLineArg("Ncell",4);
   makeShape(shp,Nc,Nc,Nc,lo,hi,comLineArg("shape",string("MPM")));
   pch=new myPatch(tmi,cst,shp);
}


