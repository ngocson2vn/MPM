// 2011 Philip Wallstedt

#include <fstream>
#include "main.h"
#include "constit.h"
#include "io.h"
using namespace std;

struct myPatch:public patch{
   ofstream hf;
   const shapeSC&shp;
   bool afterStep();
   myPatch(timeIntPtr&,constitPtr&cst,shapePtr&s);
};

bool myPatch::afterStep(){
   if(cartRank==0)for(indexT i=0;i<Nnode();++i){
      //if(massN[i]>.0001)hf<<shp.gridPosN(i)<<tab<<mag(velN[i]*massN[i])<<nwl;
      hf<<shp.gridPosN(i)<<tab<<mag(velN[i]*massN[i])<<nwl;
      //hf<<shp.gridPosN(i)<<tab<<(shp.neighborMask[i]!=0?1:0)<<nwl;
      //cerr<<i<<tab<<massN[i]<<nwl;
   }
   report.param("n0",massN[0]);
   report.param("n57",massN[57]);
   report.param("n58",massN[58]);
   report.param("n65",massN[65]);
   report.param("n114",massN[114]);
   return false;
}

myPatch::myPatch(timeIntPtr&tmi,constitPtr&cst,shapePtr&svar):shp(*svar){
   fillIndic(*this,*svar,shp.regionLo,shp.regionHi,boxIndic(shp.regionLo,shp.regionHi),vec3(2,2,2),false);
   cst=new neoHookean(1000.,.25,1.,velGradP,defGradP,volumeP,refVolP,stressP);
   makeTimeInt(tmi,this,cst,svar,comLineArg("tInt",string("CD")),0.,1.);
   for(indexT i=0;i<Npart();i+=1){
      massP[i]=.125;
      velP[i]=unit(vec3(1,1,1));
   }
   if(cartRank==0)hf.open("history.xls");  hf.precision(14);
}

void initRun(patchPtr&pch,constitPtr&cst,timeIntPtr&tmi,shapePtr&shp){
   const vec3 lo(0.,0.,0.);
   const vec3 hi(1.,1.,1.);
   int Nc=comLineArg("Ncell",8);
   makeShape(shp,Nc,Nc,Nc,lo,hi,comLineArg("shape",string("GIMP")));
   pch=new myPatch(tmi,cst,shp);
}


