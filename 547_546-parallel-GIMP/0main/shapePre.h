// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef SHAPE_PRE_H
#define SHAPE_PRE_H

#include"managedArray.h"

class arrayBuffer{
   voidPtr root;
   indexT Nbyte;
   indexT mark;
public:
   voidPtr addr;
   int want;
   arrayBuffer(){root=NULL;Nbyte=0;mark=0;}
   ~arrayBuffer(){if(root!=NULL)free(root);}
   void*operator&(){return root;}
   void checkSize(const indexT wantByteSize){
      if(Nbyte<wantByteSize){
         //cerr<<cartRank<<" enlarging from "<<Nbyte<<" to "<<wantByteSize<<endl;
         if(root!=NULL)free(root);
         root=malloc(wantByteSize);
         assert(root!=NULL);
         Nbyte=wantByteSize;
         addr=root;
         report.param("MPI buffer bytes",Nbyte);
      }
      //else cerr<<"have "<<Nbyte<<" need "<<wantByteSize<<endl;
   }
   template<typename T>void put(const T&dat){
      assert(mark*sizeof(T) < Nbyte);
      T*datptr=static_cast<T*>(root);
      datptr+=mark;
      ++mark;
      *datptr=dat;
   }
   template<typename T>T get(){
      assert(mark*sizeof(T) < Nbyte);
      T*datptr=static_cast<T*>(root);
      datptr+=mark;
      ++mark;
      return *datptr;
   }
   indexT Nput(){return mark;}
   indexT bytes(){return Nbyte;}
   void reset(){mark=0;addr=root;want=0;}
   void setbuff(const realT x){
      realT*u=static_cast<realT*>(root);
      for(indexT i=0;i<Nbyte/sizeof(realT);++i)u[i]=x;
   }
   void dumpbuff(){
      for(int n=0;n<cartSize;++n){
         //MPI_Barrier(cartComm);
         if(cartRank==n){
            realT*u=static_cast<realT*>(root);
            cerr<<cartRank<<" buffer dump "<<Nbyte<<" ==> ";
            for(indexT i=0;i<Nbyte/sizeof(realT);++i)cerr<<u[i]<<' ';
            cerr<<endl;
         }
         //MPI_Barrier(cartComm);
      }
   }
};

struct buffer{
   arrayBuffer send,recv;
   MPI_Request Sreq,Rreq;
   MPI_Status Sstat,Rstat;
   const int nbrRank;
   buffer(int d):nbrRank(d){}
};

class shapePre{
   nodeArray<vec3>curPosN;
public:
   const vec3 regionLo,regionHi,cell;
   const indexT Icol,Jrow,Klay;
   const vec3&gridPosN(const indexT i)const{return curPosN[i];}
   bool inRegion(const vec3&p)const{return regionLo<=p&&p<regionHi;}
   void setBufferSizes();
   realT dot(const nodeArray<vec3>&a,const nodeArray<vec3>&b)const;
   nodeArray<indexT>neighborMask; // should be protected
protected:
   struct contrib{
      vec3 G;
      realT w;
      indexT p,i;
      contrib(const indexT P,const indexT I,const realT W,const vec3 GG){p=P;i=I;w=W;G=GG;}
   };
   vector<contrib>conInner,conOuter;

   mutable map<indexT,buffer>allBuffers; // <neighbor mask , buffer>
   //nodeArray<indexT>neighborMask,contribCount;
   nodeArray<indexT>contribCount;
   shapePre(const indexT Nx,const indexT Ny,const indexT Nz,const vec3 lo,const vec3 hi);

   // These functions find the cell within which a particle is contained
   // The lower left node of the cell is always defined as the "parent" node of that cell.
   realT rsx(const realT x)const{return(x-regionLo.x)/cell.x+1;}
   realT rsy(const realT y)const{return(y-regionLo.y)/cell.y+1;}
   realT rsz(const realT z)const{return(z-regionLo.z)/cell.z+1;}
   indexT rsi(const realT x)const{return indexT(floor(rsx(x)));}
   indexT rsj(const realT y)const{return indexT(floor(rsy(y)));}
   indexT rsk(const realT z)const{return indexT(floor(rsz(z)));}
   indexT cellIndices(const vec3&p,indexT&i,indexT&j,indexT&k)const{ // find the cell into which x,y,z falls
      i=rsi(p.x);
      j=rsj(p.y);
      k=rsk(p.z);
      indexT nbr=0;
      if(i<1     )nbr |= 32;
      if(i>Icol-3)nbr |= 16;
      if(j<1     )nbr |= 8;
      if(j>Jrow-3)nbr |= 4;
      if(k<1     )nbr |= 2;
      if(k>Klay-3)nbr |= 1;
      //n=k*Icol*Jrow+j*Icol+i;
      return nbr;
   }

   template<typename T>void sendOuter(nodeArray<T>&gu)const{
      map<indexT,buffer>::iterator it;
      for(it=allBuffers.begin();it!=allBuffers.end();++it)it->second.send.reset();
      for(indexT n=0;n<Nnode();++n){
         const indexT self=neighborMask[n];
         if(self==0)continue; // skip interior nodes
         for(it=allBuffers.begin();it!=allBuffers.end();++it){ // check all available buffers
            const indexT proc=it->first;
            if(proc == (proc & self)){
               it->second.send.put(gu[n]);
            }
         }
      }
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         buffer&bf=it->second;
         //cerr<<"send recv "<<bf.send.Nput()<<tab<<bf.send.Nput()*sizeof(realT)<<tab<<bf.send.bytes()<<tab<<bf.recv.bytes()<<endl;
         MPI_Irecv(&bf.recv,bf.send.Nput(),*MPItype(T()),bf.nbrRank,bf.nbrRank,cartComm,&bf.Rreq);
         MPI_Isend(&bf.send,bf.send.Nput(),*MPItype(T()),bf.nbrRank,cartRank  ,cartComm,&bf.Sreq);
      }
   }

   template<typename T>void recvAndSum(nodeArray<T>&gu)const{
      map<indexT,buffer>::iterator it;
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         buffer&bf=it->second;
         waitCheckMPI(&bf.Sreq,&bf.Sstat);
         waitCheckMPI(&bf.Rreq,&bf.Rstat);
         bf.recv.reset();
      }
      for(indexT n=0;n<Nnode();++n){
         T val=gu[n];
         const indexT self=neighborMask[n];
         if(self==0)continue; // skip interior nodes
         for(it=allBuffers.begin();it!=allBuffers.end();++it){ // check all available buffers
            const indexT proc=it->first;
            if(proc == (proc & self)){
               val+=it->second.recv.get<T>(); // sum from each buffer that has a contribution
            }
         }
         gu[n]=val;
      }
   }

public:
   void sendRecvParts(const partArray<vec3>&curPosP){
      int partTypeSize;
      MPI_Type_size(MPIparticle,&partTypeSize);
      map<indexT,buffer>::iterator it;

      // sanitize buffers
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         it->second.send.reset();
         it->second.recv.reset();
      }

      // count particles to be sent to each neighbor
      for(indexT p=0;p<Npart();p++){
         indexT i,j,k;
         const indexT nbr=cellIndices(curPosP[p],i,j,k); // not storing this; should we?
         if(nbr!=0){                 // is not in local patch
            it=allBuffers.find(nbr); // is there a neighbor patch associated with the mask?
            if(it==allBuffers.end())throw stop("sendRecvParts 1 outside findable region"); // no neighbor exists (must be outside domain)
            ++it->second.send.want;  // make room for the particle to be transferred
         }
      }

      // communicate transfer counts
      for(it=allBuffers.begin();it!=allBuffers.end();++it){ // exchange expected buffer sizes
         buffer&bf=it->second;
         MPI_Irecv(&bf.recv.want,1,MPI_UNSIGNED,bf.nbrRank,bf.nbrRank,cartComm,&bf.Rreq);
         MPI_Isend(&bf.send.want,1,MPI_UNSIGNED,bf.nbrRank,cartRank  ,cartComm,&bf.Sreq);
      }

      // set send buffer sizes while waiting for recv buffer sizes
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         it->second.send.checkSize(it->second.send.want*partTypeSize);
      }

      // complete buffer size transfers
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         buffer&bf=it->second;
         waitCheckMPI(&bf.Sreq,&bf.Sstat);
         waitCheckMPI(&bf.Rreq,&bf.Rstat);
      }

      // set recv buffer sizes
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         it->second.recv.checkSize(it->second.recv.want*partTypeSize);
      }

      // copy sendable parts to buffers, and simultaneously delete them from partArrays
      indexT p=0;
      while(p<Npart()){
         indexT i,j,k;
         const indexT nbr=cellIndices(curPosP[p],i,j,k); // should be same as above
         if(nbr!=0){
            it=allBuffers.find(nbr);
            if(it==allBuffers.end())throw stop("sendRecvParts 2 outside findable region");
            globalPartArrays.partToBuff(p,it->second.send.addr);
         }
         else ++p;
      }

      // transmit parts to/from neighbors
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         buffer&bf=it->second;
         MPI_Irecv(&bf.recv,bf.recv.want,MPIparticle,bf.nbrRank,bf.nbrRank,cartComm,&bf.Rreq);
         MPI_Isend(&bf.send,bf.send.want,MPIparticle,bf.nbrRank,cartRank  ,cartComm,&bf.Sreq);
      }

      // complete particle transfers
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         buffer&bf=it->second;
         waitCheckMPI(&bf.Sreq,&bf.Sstat); // could do this after particle unpacking
         waitCheckMPI(&bf.Rreq,&bf.Rstat);
      }

      // unpack buffers to part Arrays
      for(it=allBuffers.begin();it!=allBuffers.end();++it){
         buffer&bf=it->second;
         for(int i=0;i<bf.recv.want;++i){
            globalPartArrays.partFromBuff(it->second.recv.addr);
            //cerr<<"proc "<<cartRank<<" receives from "<<it->second.nbrRank<<endl;
         }
      }
   }
};

#endif
