// Philip Wallstedt 2004-2009

#include <cstdlib> // for rand
#include <cmath> // for fmod,sin,cos
#include <cassert>
#include "func.h"
using namespace std;

// Utility for processing command line arguments
argMap::argMap(int argc,char**argv){
   for(int arg=1;arg<argc;arg+=1){
      string argstr(argv[arg]);
      string tag=argstr.substr(0,argstr.find('='));    // exception if = not found
      string val(argstr.substr(tag.size()+1));     // not sure why this can't be const
      if(find(tag)==end())insert(make_pair(tag,val));
      else throw"duplicate command line tag encountered";
   }
}

// add a particle
inline void addPart(patch&pch,const Vector2&pt,const Vector2&ps,const int Ns){
   const double vol=pch.thick*ps.x*ps.y;
   const int p=pch.appendPart();
   pch.px[p]=pt;
   pch.pX[p]=pt;
   pch.pV[p]=vol;
   pch.pm[p]=vol*pch.dens;
   pch.pF[p]=I2();
   pch.pv[p]=0.;
   pch.pCon[p].portionArray.resize(Ns);
}

// This is rarely used and not needed now
void fillRefinedSurface(patch&pch,bool(*in)(const Vector2&pos),const int Ippe,const int Sppe,const int Ns){
   assert(Sppe>=Ippe);
   const Vector2 Ips(pch.dx/double(Ippe),pch.dy/double(Ippe));
   const Vector2 Sps(pch.dx/double(Sppe),pch.dy/double(Sppe));
   for(int n=0;n<pch.Nnode();++n){
      if(!pch.inRegion(pch.gx[n]+Ips))continue;
      if(in(pch.gx[n])&&in(pch.gx[n+1])&&in(pch.gx[n+pch.I+1])&&in(pch.gx[n+pch.I])){
         for(int i=0;i<Ippe;++i){
            for(int j=0;j<Ippe;++j){
               const Vector2 ns(double(i)+.5,double(j)+.5);
               const Vector2 pt(pch.gx[n]+Ips*ns);
               addPart(pch,pt,Ips,Ns);
            }
         }
      }
      else{
         for(int i=0;i<Sppe;++i){
            for(int j=0;j<Sppe;++j){
               const Vector2 ns(double(i)+.5,double(j)+.5);
               const Vector2 pt(pch.gx[n]+Sps*ns);
               if(in(pt))addPart(pch,pt,Sps,Ns);
            }
         }
      }
   }
}

void fillAnnulusRadial(patch&pch,const Vector2&o,const double ri,const double ro,const int Ns){
   assert(pch.dx==pch.dy);
   const double dr0=pch.dx/pch.ppe;             // estimate of radius increment
   const int nr=int(Pceil((ro-ri)/dr0));         // number of rings
   const double dr=(ro-ri)/double(nr);          // true radius increment
   for(int r=0;r<nr;++r){                       // for each ring
      const double R=ri+dr*(double(r)+.5);      // R is the centerline of each ring
      const double qc=.5*pi*R;                  // circumference of quarter circle
      const int nc=4*int(Pceil(qc/dr0));         // number of segments for full circle
      const double dc=qc/double(nc/4);          // true segment circ
      const double da=2.*pi/(double(nc));       // true angle increment
      for(int c=0;c<nc;++c){                    // for each segment within a ring
         const double a=da*(double(c)+.5);      // angle of segment
         const Vector2 pt=o+Vector2(R*cos(a),R*sin(a)); // x-y coords of segment position
         if(!pch.inRegion(pt))continue;         // reject if not in bounds
         const double vol=dc*dr;                // segment volume
         const int p=pch.appendPart();
         pch.px[p]=pt;
         pch.pX[p]=pt;
         pch.pV[p]=vol;
         pch.pm[p]=vol*pch.dens;
         pch.pF[p]=I2();
         pch.pv[p]=0.;
         pch.pCon[p].portionArray.resize(Ns);
      }
   }
}

void fillAnnulusRegular(patch&pch,const Vector2&o,const double ri,const double ro,const double ppe,const int Ns){
// particles are regularly spaced within a rectangular region, whether
// they line up with grid boundaries or not.
   assert(pch.dx==pch.dy);
   const int nx=int(Pceil(2.*ro/(pch.dx/ppe)));
   const int ny=int(Pceil(2.*ro/(pch.dy/ppe)));
   const Vector2 ps(2.*ro/double(nx),2.*ro/double(ny));
   const double vol=pch.thick*ps.x*ps.y;
   for(int j=0;j<ny;j+=1){
      for(int i=0;i<nx;i+=1){
         const Vector2 ns(double(i)+.5,double(j)+.5);
         const Vector2 pt=(o-ro)+ps*ns;
         if(!pch.inRegion(pt))continue;
         if(radius(o,pt)<ri           ||radius(o,pt)>ro)continue;
         const int p=pch.appendPart();
         pch.px[p]=pt;
         pch.pX[p]=pt;
         pch.pV[p]=vol;
         pch.pm[p]=vol*pch.dens;
         pch.pF[p]=I2();
         pch.pv[p]=0.;
         pch.pCon[p].portionArray.resize(Ns);
      }
   }
}

void fillRectangleRegular(patch&pch,const Vector2&b,const Vector2&e,const double ppe,const int Ns){
// particles are regularly spaced within a rectangular region, whether
// they line up with grid boundaries or not.
   const int nx=int(Pceil((e.x-b.x)/(pch.dx/ppe)));
   const int ny=int(Pceil((e.y-b.y)/(pch.dy/ppe)));
   const Vector2 ps((e.x-b.x)/double(nx),(e.y-b.y)/double(ny));
   const double vol=pch.thick*ps.x*ps.y;
   for(int j=0;j<ny;j+=1){
      for(int i=0;i<nx;i+=1){
         const Vector2 ns(double(i)+.5,double(j)+.5);
         const Vector2 pt=b+ps*ns;
         if(!pch.inRegion(pt))continue;
         const int p=pch.appendPart();
         pch.px[p]=pt;
         pch.pX[p]=pt;
         pch.pV[p]=vol;
         pch.pm[p]=vol*pch.dens;
         pch.pF[p]=I2();
         pch.pv[p]=0.;
         pch.pCon[p].portionArray.resize(Ns);
      }
   }
}















