// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#include "main.h"
using namespace std;

class MPM:public shapeSC{
   realT S(const realT&x,const realT&h){return 1.-abs(x)/h;}
   realT G(const realT&x,const realT&h){return   -sgn(x)/h;}
   void setNeighMask(){
      neighborMask.assign(neighborMask.size(),0);
      for(indexT k=1;k<Klay-1;++k){
         for(indexT j=1;j<Jrow-1;++j){
            for(indexT i=1;i<Icol-1;++i){
               indexT nbr=0;
               if(i==1     )nbr |= 32;
               if(i==Icol-2)nbr |= 16;
               if(j==1     )nbr |= 8;
               if(j==Jrow-2)nbr |= 4;
               if(k==1     )nbr |= 2;
               if(k==Klay-2)nbr |= 1;
               neighborMask[k*Icol*Jrow+j*Icol+i]=nbr;
            }
         }
      }
   }
   void setCon(const vec3&curPos,const indexT p,const indexT n){
      for(indexT k=0;k<2;++k){
         for(indexT j=0;j<2;++j){
            for(indexT i=0;i<2;++i){
               const indexT m=n+k*Icol*Jrow+j*Icol+i;
               const vec3 r(curPos-gridPosN(m));
               const realT Sx=S(r.x,cell.x);
               const realT Sy=S(r.y,cell.y);
               const realT Sz=S(r.z,cell.z);
               const realT Gx=G(r.x,cell.x);
               const realT Gy=G(r.y,cell.y);
               const realT Gz=G(r.z,cell.z);
               const realT w=Sx*Sy*Sz;
               const realT x=Gx*Sy*Sz;
               const realT y=Gy*Sx*Sz;
               const realT z=Gz*Sx*Sy;
               if(neighborMask[m]==0)conInner.push_back(contrib(p,m,w,vec3(x,y,z)));
               else                  conOuter.push_back(contrib(p,m,w,vec3(x,y,z)));
            }
         }
      }
   }
public:
   MPM(const indexT Nx,const indexT Ny,const indexT Nz,const vec3 lo,const vec3 hi):shapeSC(Nx,Ny,Nz,lo,hi){setNeighMask();}
   indexT inCell(const vec3&p,indexT&n)const{
      indexT i,j,k;
      const indexT nbr=cellIndices(p,i,j,k);
      n=k*Icol*Jrow+j*Icol+i;
      return nbr;
   }
   void updateContribList(const patch&pch){
      conOuter.clear();
      conInner.clear();
      for(indexT p=0;p<Npart();p++){
         indexT n;
         const indexT nbr=inCell(pch.curPosP[p],n);
         if(nbr==0)setCon(pch.curPosP[p],p,n);
         else throw stop("outside findable region");
      }
      report.param("outer connectivity table size",conOuter.size());
      report.param("inner connectivity table size",conInner.size());
   }
};

class GIMP:public shapeSC{
   realT S(const realT&x,const realT&h,const realT&l){
      const realT r=abs(x);
      if(r<l)  {return 1.-(r*r+l*l)/(2.*h*l);}
      if(r<h-l){return 1.-r/h;}
      if(r<h+l){return (h+l-r)*(h+l-r)/(4.*h*l);}
      return 0.;
   }
   realT G(const realT&x,const realT&h,const realT&l){
      const realT r=abs(x);
      if(r<l  ){return -x/(h*l);}
      if(r<h-l){return -sgn(x)/h;}
      if(r<h+l){return (h+l-r)/(-2.*sgn(x)*h*l);}
      return 0.;
   }
   void setNeighMask(){
      const indexT nlay=3;
      neighborMask.assign(neighborMask.size(),0);
      for(indexT k=0;k<Klay;++k){
         for(indexT j=0;j<Jrow;++j){
            for(indexT i=0;i<Icol;++i){
               indexT nbr=0;
               if(i <  nlay     )nbr |= 32;
               if(i >= Icol-nlay)nbr |= 16;
               if(j <  nlay     )nbr |= 8;
               if(j >= Jrow-nlay)nbr |= 4;
               if(k <  nlay     )nbr |= 2;
               if(k >= Klay-nlay)nbr |= 1;
               neighborMask[k*Icol*Jrow+j*Icol+i]=nbr;
            }
         }
      }
   }
   void setCon(const vec3&curPos,const vec3&halfLenP,const indexT p,const indexT n){
      for(indexT k=0;k<3;++k){
         for(indexT j=0;j<3;++j){
            for(indexT i=0;i<3;++i){
               const indexT m=n+k*Icol*Jrow+j*Icol+i;
               const vec3 r(curPos-gridPosN(m));
               const realT Sx=S(r.x,cell.x,halfLenP.x);
               const realT Sy=S(r.y,cell.y,halfLenP.y);
               const realT Sz=S(r.z,cell.z,halfLenP.z);
               const realT Gx=G(r.x,cell.x,halfLenP.x);
               const realT Gy=G(r.y,cell.y,halfLenP.y);
               const realT Gz=G(r.z,cell.z,halfLenP.z);
               const realT w=Sx*Sy*Sz;
               const realT x=Gx*Sy*Sz;
               const realT y=Gy*Sx*Sz;
               const realT z=Gz*Sx*Sy;
               if(neighborMask[m]==0)conInner.push_back(contrib(p,m,w,vec3(x,y,z)));
               else                  conOuter.push_back(contrib(p,m,w,vec3(x,y,z)));
            }
         }
      }
   }
public:
   GIMP(const indexT Nx,const indexT Ny,const indexT Nz,const vec3 lo,const vec3 hi):shapeSC(Nx,Ny,Nz,lo,hi){setNeighMask();}
   indexT inCell(const vec3&p,indexT&n)const{
      indexT i,j,k;
      const indexT nbr=cellIndices(p,i,j,k);
      if(nbr!=0)return nbr;
      const realT enx=rsx(p.x)-realT(i);
      const realT eny=rsy(p.y)-realT(j);
      const realT enz=rsz(p.z)-realT(k);
      const indexT io=(enx<.5?i-1:i);
      const indexT jo=(eny<.5?j-1:j);
      const indexT ko=(enz<.5?k-1:k);
      n=ko*Icol*Jrow+jo*Icol+io;
      return nbr;
   }
   void updateContribList(const patch&pch){
      conOuter.clear();
      conInner.clear();
      for(indexT p=0;p<Npart();p++){
         indexT n;
         const indexT nbr=inCell(pch.curPosP[p],n);
         if(nbr==0)setCon(pch.curPosP[p],pch.halfLenP[p],p,n);
         else throw stop("outside findable region");
      }
      report.param("outer connectivity table size",conOuter.size());
      report.param("inner connectivity table size",conInner.size());
   }
};

class spline:public shapeSC{
   realT S(const realT&x,const realT&h){
      const realT r=abs(x);
      if(r< .5*h){return .75-r*r/(h*h);}
      if(r<1.5*h){return (1.5*h-r)*(1.5*h-r)/(2.*h*h);}
      return 0.;
   }
   realT G(const realT&x,const realT&h){
      const realT r=abs(x);
      if(r< .5*h){return -x/(.5*h*h);}
      if(r<1.5*h){return (1.5*h-r)/(-sgn(x)*h*h);}
      return 0.;
   }
   void setNeighMask(){
      const indexT nlay=3;
      neighborMask.assign(neighborMask.size(),0);
      for(indexT k=0;k<Klay;++k){
         for(indexT j=0;j<Jrow;++j){
            for(indexT i=0;i<Icol;++i){
               indexT nbr=0;
               if(i <  nlay     )nbr |= 32;
               if(i >= Icol-nlay)nbr |= 16;
               if(j <  nlay     )nbr |= 8;
               if(j >= Jrow-nlay)nbr |= 4;
               if(k <  nlay     )nbr |= 2;
               if(k >= Klay-nlay)nbr |= 1;
               neighborMask[k*Icol*Jrow+j*Icol+i]=nbr;
            }
         }
      }
   }
   void setCon(const vec3&curPos,const indexT p,const indexT n){
      for(indexT k=0;k<3;++k){
         for(indexT j=0;j<3;++j){
            for(indexT i=0;i<3;++i){
               const indexT m=n+k*Icol*Jrow+j*Icol+i;
               const vec3 r(curPos-gridPosN(m));
               const realT Sx=S(r.x,cell.x);
               const realT Sy=S(r.y,cell.y);
               const realT Sz=S(r.z,cell.z);
               const realT Gx=G(r.x,cell.x);
               const realT Gy=G(r.y,cell.y);
               const realT Gz=G(r.z,cell.z);
               const realT w=Sx*Sy*Sz;
               const realT x=Gx*Sy*Sz;
               const realT y=Gy*Sx*Sz;
               const realT z=Gz*Sx*Sy;
               if(neighborMask[m]==0)conInner.push_back(contrib(p,m,w,vec3(x,y,z)));
               else                  conOuter.push_back(contrib(p,m,w,vec3(x,y,z)));
            }
         }
      }
   }
public:
   spline(const indexT Nx,const indexT Ny,const indexT Nz,const vec3 lo,const vec3 hi):shapeSC(Nx,Ny,Nz,lo,hi){setNeighMask();}
   indexT inCell(const vec3&p,indexT&n)const{
      indexT i,j,k;
      const indexT nbr=cellIndices(p,i,j,k);
      if(nbr!=0)return nbr;
      const realT enx=rsx(p.x)-realT(i);
      const realT eny=rsy(p.y)-realT(j);
      const realT enz=rsz(p.z)-realT(k);
      const indexT io=(enx<.5?i-1:i);
      const indexT jo=(eny<.5?j-1:j);
      const indexT ko=(enz<.5?k-1:k);
      n=ko*Icol*Jrow+jo*Icol+io;
      return nbr;
   }
   void updateContribList(const patch&pch){
      conOuter.clear();
      conInner.clear();
      for(indexT p=0;p<Npart();p++){
         indexT n;
         const indexT nbr=inCell(pch.curPosP[p],n);
         if(nbr==0)setCon(pch.curPosP[p],p,n);
         else throw stop("outside findable region");
      }
      report.param("outer connectivity table size",conOuter.size());
      report.param("inner connectivity table size",conInner.size());
   }
};

void makeShape(shapePtr&shp,const indexT Nx,const indexT Ny,const indexT Nz,const vec3 lo,const vec3 hi,string s){
   int dims[3]={0,0,0};
   int periods[3]={0,0,0};
   int coords[3]={0,0,0};
   MPI_Cart_get(cartComm,3,dims,periods,coords);
   const vec3 nproc(dims[0],dims[1],dims[2]);
   const vec3 myCoords(coords[0],coords[1],coords[2]);
   const vec3 Ncell=vec3(Nx,Ny,Nz);
   const vec3 cell((hi-lo)/Ncell);
   const vec3 loIdx=round( myCoords    *Ncell/nproc);
   const vec3 hiIdx=round((myCoords+1.)*Ncell/nproc);
   const vec3 NmyCell=hiIdx-loIdx;
   const vec3 myLo=lo+cell*loIdx;
   const vec3 myHi=lo+cell*hiIdx;
   const indexT myNx=indexT(NmyCell.x);
   const indexT myNy=indexT(NmyCell.y);
   const indexT myNz=indexT(NmyCell.z);
   shp=NULL;
   if(s=="MPM"   ){shp=new MPM   (myNx,myNy,myNz,myLo,myHi);}
   if(s=="GIMP"  ){shp=new GIMP  (myNx,myNy,myNz,myLo,myHi);}
   if(s=="spline"){shp=new spline(myNx,myNy,myNz,myLo,myHi);}
   if(shp==NULL)throw stop("bad shape string");
   shp->setBufferSizes();


}





