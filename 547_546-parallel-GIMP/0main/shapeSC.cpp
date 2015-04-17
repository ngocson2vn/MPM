// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include"main.h"
using namespace std;

realT shapePre::dot(const nodeArray<vec3>&a,const nodeArray<vec3>&b)const{
   const realT*aptr=&a[0].x;
   const realT*bptr=&b[0].x;
   complex<realT>sum(0.,0.);
   complex<realT>addend(0.,0.);
   for(indexT j=0;j<3*Nnode();++j){
      addend.real()=(aptr[j]*bptr[j])/realT(contribCount[j/3]);
      DDadd(&addend.real(),&sum.real(),NULL,NULL);
   }
   MPI_Allreduce(MPI_IN_PLACE,&sum,1,MPI::DOUBLE_COMPLEX,MPI_DDadd,cartComm);
   assert(sum.real() == sum.real());
   return sum.real();
}

shapePre::shapePre(const indexT Nx,const indexT Ny,const indexT Nz,const vec3 lo,const vec3 hi):
               regionLo(lo),
               regionHi(hi),
               cell((hi-lo)/vec3(realT(Nx),realT(Ny),realT(Nz))),
               Icol(Nx+3), // must keep full boundary layer for GIMP
               Jrow(Ny+3),
               Klay(Nz+3){
   globalNodeArrays.resizeAll(Icol*Jrow*Klay);
   for(int k=0;k<int(Klay);++k){
      const realT z=(k-1)*cell.z+regionLo.z;
      for(int j=0;j<int(Jrow);++j){
         const realT y=(j-1)*cell.y+regionLo.y;
         for(int i=0;i<int(Icol);++i){
            const realT x=(i-1)*cell.x+regionLo.x;
            curPosN[k*Icol*Jrow+j*Icol+i]=vec3(x,y,z);
         }
      }
   }
   report.param("Nnode",Nnode());
   report.param("regionLo",regionLo);
   report.param("regionHi",regionHi);
   report.param("cell sizes",cell);
   report.param("Icols",Icol);
   report.param("Jrows",Jrow);
   report.param("Klayers",Klay);

   int dims[3],periods[3],coords[3];
   MPI_Cart_get(cartComm,3,dims,periods,coords);
   const int I=coords[0];
   const int J=coords[1];
   const int K=coords[2];
   const int Iproc=dims[0];
   const int Jproc=dims[1];
   const int Kproc=dims[2];
   // my coords are I J K
   // total number of procs in domain is Iproc Jproc Kproc
   for(int k=K-1;k<K+2;++k)if(k>=0 && k<Kproc){
      for(int j=J-1;j<J+2;++j)if(j>=0 && j<Jproc){
         for(int i=I-1;i<I+2;++i)if(i>=0 && i<Iproc){
            int nbrRank;
            int nbrCoords[]={i,j,k};
            MPI_Cart_rank(cartComm,nbrCoords,&nbrRank);
            if(nbrRank!=cartRank){
               indexT nbr=0;
               if(i<I)nbr |= 32;
               if(i>I)nbr |= 16;
               if(j<J)nbr |= 8;
               if(j>J)nbr |= 4;
               if(k<K)nbr |= 2;
               if(k>K)nbr |= 1;
               assert(nbr!=0);
               allBuffers.insert(make_pair(nbr,buffer(nbrRank)));
            }
         }
      }
   }
}

void shapePre::setBufferSizes(){
   map<indexT,indexT>bcount; // <up to 26 buffers , N nodes in the buffer>
   for(map<indexT,buffer>::iterator it=allBuffers.begin();it!=allBuffers.end();++it)bcount[it->first]=0;
   for(indexT n=0;n<Nnode();++n){
      const indexT self=neighborMask[n];
      contribCount[n]=1;     // at least one patch owns me
      if(self==0)continue;   // skip interior nodes
      for(map<indexT,indexT>::iterator it=bcount.begin();it!=bcount.end();++it){ // check all available buffers
         const indexT proc=it->first;
         if(proc == (proc & self)){
            it->second+=1;   // count nodes in buffer
            ++contribCount[n]; // count times node has contributed
         }
      }
      if(contribCount[n]==1){
         neighborMask[n]=0; // removing nodes on domain surface
      }
   }
   for(map<indexT,buffer>::iterator it=allBuffers.begin();it!=allBuffers.end();++it){
      map<indexT,indexT>::iterator bc=bcount.find(it->first);
      assert(bc!=bcount.end());
      indexT bsize=bc->second;
      it->second.send.checkSize(bsize*indexT(sizeof(vec3)));
      it->second.recv.checkSize(bsize*indexT(sizeof(vec3)));
   }
}

void shapeSC::integrate(const partArray<realT>&pu,nodeArray<realT>&gu)const{
   gu.assign(gu.size(),0.);
   for(indexT i=0;i<conOuter.size();++i){
      const contrib&cn=conOuter[i];
      gu[cn.i]+=pu[cn.p]*cn.w;
   }
   sendOuter(gu);
   for(indexT i=0;i<conInner.size();++i){
      const contrib&cn=conInner[i];
      gu[cn.i]+=pu[cn.p]*cn.w;
   }
   recvAndSum(gu);
}

void shapeSC::integrate(const partArray<vec3>&pu,nodeArray<vec3>&gu)const{
   gu.assign(gu.size(),vec3(0.,0.,0.));
   for(indexT i=0;i<conOuter.size();++i){
      const contrib&cn=conOuter[i];
      gu[cn.i]+=pu[cn.p]*cn.w;
   }
   sendOuter(gu);
   for(indexT i=0;i<conInner.size();++i){
      const contrib&cn=conInner[i];
      gu[cn.i]+=pu[cn.p]*cn.w;
   }
   recvAndSum(gu);
}

void shapeSC::divergence(const partArray<mat3>&stressP,const partArray<realT>&volumeP,nodeArray<vec3>&forceN)const{
   forceN.assign(forceN.size(),vec3(0.,0.,0.));
   for(indexT i=0;i<conOuter.size();++i){
      const contrib&cn=conOuter[i];
      forceN[cn.i]-=(volumeP[cn.p]*stressP[cn.p]).inner(cn.G);
   }
   sendOuter(forceN);
   for(indexT i=0;i<conInner.size();++i){
      const contrib&cn=conInner[i];
      forceN[cn.i]-=(volumeP[cn.p]*stressP[cn.p]).inner(cn.G);
   }
   recvAndSum(forceN);
}

void shapeSC::interpolate(partArray<vec3>&pu,const nodeArray<vec3>&gu)const{
   pu.assign(pu.size(),vec3(0.,0.,0.));
   for(indexT i=0;i<conOuter.size();++i){
      const contrib&cn=conOuter[i];
      pu[cn.p]+=gu[cn.i]*cn.w;
   }
   for(indexT i=0;i<conInner.size();++i){
      const contrib&cn=conInner[i];
      pu[cn.p]+=gu[cn.i]*cn.w;
   }
}

void shapeSC::gradient(partArray<mat3>&pu,const nodeArray<vec3>&gu)const{
   pu.assign(pu.size(),mat3(0.,0.,0.,0.,0.,0.,0.,0.,0.));
   for(indexT i=0;i<conOuter.size();++i){
      const contrib&cn=conOuter[i];
      pu[cn.p]+=gu[cn.i].outer(cn.G);
   }
   for(indexT i=0;i<conInner.size();++i){
      const contrib&cn=conInner[i];
      pu[cn.p]+=gu[cn.i].outer(cn.G);
   }
}

void shapeSC::gradient(partArray<vec3>&pu,const nodeArray<realT>&gu)const{
   pu.assign(pu.size(),vec3(0.,0.,0.));
   for(indexT i=0;i<conOuter.size();++i){
      const contrib&cn=conOuter[i];
      pu[cn.p]+=gu[cn.i]*cn.G;
   }
   for(indexT i=0;i<conInner.size();++i){
      const contrib&cn=conInner[i];
      pu[cn.p]+=gu[cn.i]*cn.G;
   }
}





