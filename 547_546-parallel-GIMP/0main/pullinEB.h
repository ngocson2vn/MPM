// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef PULLIN_EB_H
#define PULLIN_EB_H

#include<iostream>
#include<vector>
#include<limits>
#include<cmath>
using namespace std;

class pullinEB{
   const double beamL,beamH,beamB,Ymod,gap0,pmtv;
   double volt;
   const unsigned n;
   vector<double>pos,currDef;
public:
   pullinEB(const double a,
            const double b,
            const double c,
            const double pois,
            const double d,
            const double e,
            const double f,
            const unsigned g):
               beamL(a),
               beamH(b),
               beamB(c),
               Ymod(d/(1.-pois*pois)),
               gap0(e),
               pmtv(f),
               n(g),
               pos(n),
               currDef(n){
      for(unsigned i=0;i<n;++i)pos[i]=beamL*(.5+double(i))/double(n);
   }
   double presGivenGap(const double gap_mm,const double myvolt){ // electro-static load
      const double gap=1e-3*gap_mm;                            // (meter)
      const double fperA=pmtv*myvolt*myvolt/(2.*gap*gap);          // (Newton/meter^2)
      return 1e-6*fperA;                                       // (Newton/mm^2)
   }
   double defPoint(const double x,const double xforce,const double force){
      const double momI=beamB*beamH*beamH*beamH/12.;
      const double a=xforce;
      const double b=beamL-xforce;
      const double t1=force/(6.*Ymod*momI);
      const double t2=b*b*x*x*x*(beamL+2.*a)/(beamL*beamL*beamL);
      const double t3=3.*a*b*b*x*x/(beamL*beamL);
      const double xa=x-a;
      const double t4=(x>a?xa*xa*xa:0.);
      return t1*(t2-t3-t4);
   }
private:
   double sumForce(const double x,const double myvolt,const vector<double>&def,const vector<double>&mypos){
      double sum=0.;
      for(unsigned j=0;j<mypos.size();++j){ // sum up contributions to the point's deflection from all forces along the beam
         const double gap=gap0-def[j];
         //const double gap=gap0; // linear check
         const double fperL=beamB*presGivenGap(gap,myvolt);
         const double force=-fperL*(beamL/double(mypos.size()));
         sum+=defPoint(x,pos[j],force);
      }
      return sum;
   }
public:
   bool profileVolt(const double vl){
      volt=vl;
      vector<double>prevDef(currDef.size());
      for(unsigned i=0;i<pos.size();++i)prevDef[i]=0.;
      for(unsigned i=0;i<pos.size();++i)currDef[i]=0.;
      double norm2=gap0;
      for(unsigned k=0;k<200 && norm2>numeric_limits<double>::epsilon()*gap0;++k){
         currDef.swap(prevDef);
         norm2=0.;
         for(unsigned i=0;i<pos.size();++i){ // for each point of interest
            currDef[i]=sumForce(pos[i],volt,prevDef,pos);
            if(abs(currDef[i])>gap0){
               cerr<<"pullinEB.profileVolt: Pull-in detected!"<<endl;
               volt=0./0.;
               return false;
            }
            const double diff=currDef[i]-prevDef[i];
            norm2+=diff*diff;
         }
         //cerr<<k<<"   "<<norm2<<endl;
      }
      return true;
   }
   double defPos(const double x){return sumForce(x,volt,currDef,pos);}
};


#endif
