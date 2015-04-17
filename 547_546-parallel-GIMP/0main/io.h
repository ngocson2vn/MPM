// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef IO_H
#define IO_H

#include"main.h"

class boxIndic{
   const vec3 lo,hi;
   const bool invert;
public:
   boxIndic(const vec3&l,const vec3&h,bool iv=false):lo(l),hi(h),invert(iv){}
   bool operator()(const vec3&pt)const{
      const bool keep=pt.x<=hi.x
                   && pt.y<=hi.y
                   && pt.z<=hi.z
                   && pt.x>=lo.x
                   && pt.y>=lo.y
                   && pt.z>=lo.z;
      return (invert?!keep:keep); // Default is to keep particles located inside the box.
      // But with invert==true, we only keep particles located outside the box.
   }
};

template<typename T>
void fillIndic(patch&pch,shapeSC&shp,const vec3&b,const vec3&e,const T&indic,vec3 ppe=vec3(2,2,2),bool massNorm=false){
// Particles are regularly spaced within a box, whether they line up with grid
// boundaries or not.  The ppe is a suggestion but may not be exactly followed.
   report.param("ppd",ppe);
   const vec3 box=e-b;
   const vec3 wantSize=box/(shp.cell/ppe);
   const vec3 ppdim(round(wantSize.x),
                    round(wantSize.y),
                    round(wantSize.z));
   const int nx=int(ppdim.x);
   const int ny=int(ppdim.y);
   const int nz=int(ppdim.z);
   const vec3 partSize=box/ppdim;
   const vec3 partAreas(partSize.y*partSize.z,partSize.x*partSize.z,partSize.x*partSize.y);
   const realT vol=(box.x*box.y*box.z)/realT(nx*ny*nz);
   int count=0;
   for(int k=0;k<nz;++k){
      for(int j=0;j<ny;++j){
         for(int i=0;i<nx;++i){
            const vec3 ns(realT(i)+.5,realT(j)+.5,realT(k)+.5);
            const vec3 pt=b+box*ns/ppdim;
            if(shp.inRegion(pt) && indic(pt))++count; // count expected number of particles
         }
      }
   }
   globalPartArrays.resizeAll(count);
   count=0;
   for(int k=0;k<nz;++k){
      for(int j=0;j<ny;++j){
         for(int i=0;i<nx;++i){
            const vec3 ns(realT(i)+.5,realT(j)+.5,realT(k)+.5);
            const vec3 pt=b+box*ns/ppdim;
            if(shp.inRegion(pt) && indic(pt)){
               pch.curPosP[count]=pt;
               pch.refPosP[count]=pt;
               pch.volumeP[count]=vol;
               pch.refVolP[count]=vol;
               pch.halfLenP[count]=.5*partSize;
               pch.refHalfP[count]=.5*partSize;
               pch.massP[count]=0.;
               pch.velP[count]=0.;
               pch.defGradP[count]=I3();
               pch.surfAreaP[count]=0.;
               vec3 sn;
               // find out if neighbors would exist; if not, I am a surface particle
               if(!indic(b+box*vec3(realT(i-1)+.5,realT(j  )+.5,realT(k  )+.5)/ppdim)){sn+=vec3(-1, 0, 0);}
               if(!indic(b+box*vec3(realT(i+1)+.5,realT(j  )+.5,realT(k  )+.5)/ppdim)){sn+=vec3( 1, 0, 0);}
               if(!indic(b+box*vec3(realT(i  )+.5,realT(j-1)+.5,realT(k  )+.5)/ppdim)){sn+=vec3( 0,-1, 0);}
               if(!indic(b+box*vec3(realT(i  )+.5,realT(j+1)+.5,realT(k  )+.5)/ppdim)){sn+=vec3( 0, 1, 0);}
               if(!indic(b+box*vec3(realT(i  )+.5,realT(j  )+.5,realT(k-1)+.5)/ppdim)){sn+=vec3( 0, 0,-1);}
               if(!indic(b+box*vec3(realT(i  )+.5,realT(j  )+.5,realT(k+1)+.5)/ppdim)){sn+=vec3( 0, 0, 1);}
               if(len2(sn)>machEps){
                  pch.surfAreaP[count]=partAreas;
                  pch.surfNormP[count]=unit(sn);
               }
               else{
                  pch.surfAreaP[count]=0.;
                  pch.surfNormP[count]=0.;
               }
               ++count;
            }
         }
      }
   }
   report.param("Npart",Npart());
   if(massNorm){
      shp.updateContribList(pch);
      for(indexT i=0;i<Npart();++i)pch.massP[i]=1.;
      shp.integrate(pch.massP,pch.massN);
      for(indexT i=0;i<Npart();++i)pch.massP[i]=0.;
      shp.gradient(pch.surfNormP,pch.massN); // surfNorm based on mass gradient
      for(indexT i=0;i<Npart();++i){
         pch.surfNormP[i]=(len2(pch.surfNormP[i])>machEps ? -unit(pch.surfNormP[i]) : vec3());
      }
   }
}

inline bool partProbe(const shapeSC&shp,patch&pch,const vec3&pos){
   if(shp.inRegion(pos)){
      indexT idx=0;
      realT rmin=mag(pch.refPosP[idx]-pos);
      for(indexT i=1;i<Npart();++i){
         realT r=mag(pch.refPosP[i]-pos);
         if(r<rmin){idx=i;rmin=r;}
      }
      pch.probeP[idx]=1;
      report.param("partProbe",pch.refPosP[idx]);
      return true;
   }
   else return false;
}

inline indexT nodeProbe(const shapeSC&shp,const vec3&pos){
   indexT idx=0;
   realT rmin=mag(shp.gridPosN(idx)-pos);
   for(indexT i=1;i<Nnode();++i){
      realT r=mag(shp.gridPosN(i)-pos);
      if(r<rmin){idx=i;rmin=r;}
   }
   report.param("nodeProbe",shp.gridPosN(idx));
   return idx;
}

#endif
