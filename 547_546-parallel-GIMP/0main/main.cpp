// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include<ctime>
#include"main.h"
#include"timeInt.h"
using namespace std;

arrayGroup  globalPartArrays;
arrayGroup  globalNodeArrays;
CLargMap    globalArgMap;
paramReport report;
MPI_Comm cartComm;
int cartRank,cartSize;
MPI_Datatype MPIindexT,MPIrealT,MPIvec3,MPImat3,MPIparticle;
MPI_Op MPI_DDadd;

//void here(const int h){
   //MPI_Barrier(cartComm);
   //if(cartRank==0){
      //cerr<<"here"<<h<<' ';
      //cin.get();
   //}
   //MPI_Barrier(cartComm);
//}

void here(const int h){
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Barrier(MPI_COMM_WORLD);
   if(rank==0){
      cerr<<"here"<<h<<' ';
      cin.get();
   }
   MPI_Barrier(MPI_COMM_WORLD);
}

void crash(stop&sp){
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   report.param("crash encountered with world rank",rank);
   if(sp[0]==string("otherProcFailed"))return; // must not call allreduce
   sp.push_back(string("world rank=")+toStr(rank));
   cerr<<"Early termination requested with messages:"<<endl;
   for(indexT i=0;i<sp.size();++i)cerr<<" * "<<sp[i]<<endl;
   int ok=1;
   MPI_Allreduce(MPI_IN_PLACE,&ok,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
}

void commCheck(){
   int ok=0;
   MPI_Allreduce(MPI_IN_PLACE,&ok,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
   if(ok>0)throw stop("otherProcFailed"); // must throw, but not call allreduce again
}

void waitCheckMPI(MPI_Request*request,MPI_Status*status){
   const realT wasteLimit=1.+.1*realT(cartSize); // tenth of a second for each proc
   int flag=false;
   const clock_t startClock=clock();
   while(flag==false){
      MPI_Test(request,&flag,status);
      const realT timeWasted=realT(clock()-startClock)/CLOCKS_PER_SEC;
      if(timeWasted > wasteLimit)commCheck();
   }
}

void initMPI(){
   MPI_Type_dup(MPI_DOUBLE,&MPIrealT);
   MPI_Type_dup(MPI_UNSIGNED,&MPIindexT);
   MPI_Type_contiguous(3,MPI_DOUBLE,&MPIvec3);
   MPI_Type_commit(&MPIvec3);
   MPI_Type_contiguous(9,MPI_DOUBLE,&MPImat3);
   MPI_Type_commit(&MPImat3);
   const vec3 nproc=comLineArg("np",vec3(1,1,1));
   int dims[3]={int(nproc.x),int(nproc.y),int(nproc.z)};
   int periods[3]={0,0,0};
   MPI_Cart_create(MPI_COMM_WORLD,3,dims,periods,1,&cartComm);
   MPI_Comm_rank(cartComm,&cartRank);
   MPI_Comm_size(cartComm,&cartSize);
   MPI_Op_create(DDadd,1,&MPI_DDadd);
}

realT findMaxVel(const partArray<vec3>&v){
   realT maxVel2=0.;
   for(indexT i=0;i<Npart();++i){
      const realT v2=len2(v[i]);
      if(v2>maxVel2)maxVel2=v2; // avoiding square roots
   }
   return sqrt(maxVel2);
}

void progressIndicator(const patch&pch,const timeIntSC&tmi,const clock_t t0){
   const realT curWall=realT(clock()-t0)/CLOCKS_PER_SEC;
   const realT trat=tmi.finalTime/tmi.elapsedTime();
   report.param("max part velocity",findMaxVel(pch.velP));
   report.param("estimated wall time remaining",timeTag(curWall*(trat-1.)));
   report.param("estimated total wall time",timeTag(curWall*trat));
   report.write();
}

int main(int argc,char**argv){
   try{
      const clock_t startClock=clock();
      if(cartRank==0)cerr<<"starting . . .";
      globalArgMap.init(argc,argv);
      report.init();
      MPI_Init(NULL,NULL);
      initMPI();
      patchPtr pch=NULL;
      constitPtr cst=NULL;
      shapePtr shp=NULL;
      timeIntPtr tmi=NULL;
      initRun(pch,cst,tmi,shp);
      globalPartArrays.makeMPIbundleType(&MPIparticle);
      cst->update(0.); // time step of zero
      globalArgMap.reportUnusedTags();
      report.write();
      try{
         do{
            shp->sendRecvParts(pch->curPosP);
            shp->updateContribList(*pch);
            tmi->nextStep(*pch);
            tmi->advance(*shp);
            report.param("elapsedTime",tmi->elapsedTime());
            report.param("incCount",tmi->incCount());
            progressIndicator(*pch,*tmi,startClock);
         }while(pch->afterStep());
      }
      catch(stop&sp){sp+="Defined Exception"; crash(sp);}
      catch(const exception&error){stop sp(error.what()); sp+="Standard Exception"; crash(sp);}
      catch(...){stop sp("Unknown Exception"); crash(sp);}
      delete tmi;
      delete pch;
      delete cst;
      delete shp;
      report.param("wallTime",realT(clock()-startClock)/CLOCKS_PER_SEC);
      report.writeFinal();
      MPI_Finalize();
      if(cartRank==0)cerr<<" complete."<<endl;
   }
   catch(stop&sp){sp+="Non-loop Defined Exception"; crash(sp);}
   catch(const exception&error){stop sp(error.what()); sp+="Non-loop Standard Exception"; crash(sp);}
   catch(...){stop sp("Non-loop Unknown Exception"); crash(sp);}
}
