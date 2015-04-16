// Philip Wallstedt 2004-2009
#include <fstream>
#include <ctime>
#include "func.h"
using namespace std;

static bool single=true;

double ultimateStrength=numeric_limits<double>::max();

double planeStrainNeoHookean(const patch&pch,const Matrix2&F,Matrix2&S){
   const double v=pch.pois;
   const double l=pch.Ymod*v/((1.+v)*(1.-2.*v));
   const double m=.5*pch.Ymod/(1.+v);
   const double Ja=det(F);
   S =I2(l*log(Ja)/Ja);
   S+=m/Ja*(F.inner(trans(F))-I2());
   if(vonMises(S)>ultimateStrength)S=0.;
   return Ja;
}
double timeIntSC::getStress(const patch&pch,const Matrix2&F,Matrix2&S){return planeStrainNeoHookean(pch,F,S);}

void history(patch&pch,ostream&os,int&fc){ // colored splot
   if(single&&pch.incCount%16==0){
      fc+=1;
      os<<'#'<<fc<<'\n';

      for(int i=0;i<pch.Npart();++i){
         Matrix2 s;
         planeStrainNeoHookean(pch,pch.pF[i],s);
         //os<<pch.px[i]<<'\t'<<vonMises(s)<<'\t'<<pch.partSize/double(pch.I)<<'\n';
         os<<pch.px[i]<<'\t'<<vonMises(s)<<'\t'<<pch.partSize/double(pch.I)<<'\n';
         //os<<pch.px[i]<<'\t'<<magnitude(pch.px[i]-pch.pX[i])<<'\t'<<pch.param1/double(pch.I)<<'\n';
      }

      /* for(int i=0;i<pch.Nnode();i+=1){
         //if(pch.gx[i].x<pch.cellTol||pch.gx[i].x>1.-pch.cellTol)continue;
         //if(pch.gx[i].y<pch.cellTol||pch.gx[i].y>1.-pch.cellTol)continue;
         os<<pch.gx[i]<<'\t'<<pch.gMirr[i].y<<'\t'<<pch.param1/double(pch.I)<<'\n';
      } */

      /* for(int i=0;i<pch.Nquad();i+=1){
         if(pch.inRegion(pch.qx[i]))os<<pch.qx[i]<<'\t'<<vonMises(pch.qs[i])<<'\t'<<pch.param1/double(pch.I)<<'\n';
      } */
      os<<"\n\n";
   }
}

void timeIntSC::applyGridBC(patch&){}

int main(int argc,char**argv){
   try{
      const double CPS=double(CLOCKS_PER_SEC);
      clock_t t=clock();
      argMap am(argc,argv); // assemble command line parameter map
      int Nc=20;                       am["Ncell"]>>Nc;
      string shpstr="GIMP";            am["shape"]>>shpstr;
      string tmistr="cen";             am["tInt" ]>>tmistr;
      double CFL=.2;                   am["CFL"  ]>>CFL;
                                       am["single"]>>single;
      shapeSC&shp=makeShape(shpstr);
      patch pch(2*Nc,Nc,0.,0.,2.,1.,shp.Nghost(),1.);
      timeIntSC&ti=makeTimeInt(tmistr,shp,pch);
      pch.ppe=2.;                      am["ppe"  ]>>pch.ppe;
      pch.load=1.;                     am["load" ]>>pch.load;
      pch.load=1.;                     am["load" ]>>pch.load;
                                       am["ult"  ]>>ultimateStrength;
      pch.Ymod=1e7;
      pch.dens=1e3;
      pch.vwav=sqrt(pch.Ymod/pch.dens);
      pch.pois=.3;
      pch.damp=.0;
      pch.partSize=32.;                am["p1"   ]>>pch.partSize; // gnuplot point size
      pch.dt=min(pch.dx,pch.dy)*CFL/sqrt(pch.Ymod/pch.dens);
      fillRectangleRegular(pch,Vector2(.2+2.*pch.dx,.3),Vector2(1.+2.*pch.dx,.5),pch.ppe,shp.Nsupport()); // impacter
      fillRectangleRegular(pch,Vector2(1.+4.*pch.dx,.1),Vector2(1.2+4.*pch.dx,.9),pch.ppe,shp.Nsupport()); // target
      pch.elaps=.0/pch.vwav;
      for(int i=0;i<pch.Npart();i+=1){
         pch.px[i]=pch.pX[i];
         pch.pv[i].y=0.;
         pch.pv[i].x=(pch.px[i].x<.8+3.*pch.dx?pch.load:0.);
         if(pch.px[i].x<.8+3.*pch.dx)pch.pm[i]*=100.;
         pch.pF[i]=I2();
         pch.pfe[i]=pch.pm[i]*Vector2(0.,-9.8);
      }
      shp.updateContribList(pch);
// main loop
      ofstream hf("history.xls");
      int frameCount=0;
      double Li=0.,L2=0.,L1=0.;
      try{
         while(pch.elaps<8./pch.vwav){
         //while(pch.incCount<1){
            ti.advance(pch);
            pch.elaps+=pch.dt;
            pch.incCount+=1;
            history(pch,hf,frameCount);
            shp.updateContribList(pch);
            cerr<<'~';
         }
      }
      catch(const char*s){cerr<<"Defined Exception:"<<s<<endl;L1=L2=Li=1.;}
      catch(const std::exception&error){cerr<<"Standard Exception:"<<error.what()<<endl;}
      catch(...){cerr<<"Unknown Exception"<<endl;}
      cerr<<'\n';
// post process
      if(single){
         cerr<<"Wall time:"<<double(clock()-t)/CPS<<endl;
         cout<<frameCount;
      }
      else{
         cout<<pch.Npart()<<'\t'
             <<pch.I*pch.J<<'\t'
             <<pch.incCount<<'\t'
             <<double(clock()-t)/CPS<<'\t';
      }
   }catch(...){cerr<<"Non-loop Exception"<<endl;}
}


