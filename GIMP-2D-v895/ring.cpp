
// Philip Wallstedt 2004-2009
// This problem is taken from Wallstedt and Guilkey
// "An evaluation of explicit time integration schemes for use with the
// generalized interpolation material point method" JCP 2008.

#include <fstream>
#include <ctime>
#include "func.h"
using namespace std;

static bool single=true;
static double globalLiError=0.;
static double globalL2Error=0.;
static double globalL1Error=0.;

double zeroPois(const patch&pch,const Matrix2&F,Matrix2&S){
   const double Ja=det(F);
   S=.5*pch.Ymod/Ja*(F.inner(trans(F))-I2());
   return Ja;
}
double timeIntSC::getStress(const patch&pch,const Matrix2&F,Matrix2&S){return zeroPois(pch,F,S);}

class exact{
   int ic;
public:
   const patch&pch;
   const double c,A,E,Ri,Ro,c3,c2,c1;
   double t,S,C;
   exact(const patch&p):
         pch(p),
         c(p.vwav),
         A(p.load),
         E(p.Ymod),
         Ri(.5),
         Ro(1.),
         c3(-2.       /(Ro*Ro*(Ro-3.*Ri))),
         c2(3.*(Ri+Ro)/(Ro*Ro*(Ro-3.*Ri))),
         c1(-6.*Ri    /(   Ro*(Ro-3.*Ri))){
      ic=-1;
   }
   void base(){
      if(pch.incCount!=ic){
         ic=pch.incCount;
         t=pch.elaps;
         S=sin(c*pi*t);
         C=cos(c*pi*t);
      }
   }
   Vector2 displace(const Vector2&pos){
      const double X=pos.x;
      const double Y=pos.y;
      const double R=sqrt(X*X+Y*Y);
      base();
      const double W=A*C*(c3*R*R*R+c2*R*R+c1*R)/R;
      return Vector2(W*X,W*Y);
   }
   Vector2 velocity(const Vector2&pos){
      const double X=pos.x;
      const double Y=pos.y;
      const double R=sqrt(X*X+Y*Y);
      base();
      const double W=-pi*c*A*S*(c3*R*R*R+c2*R*R+c1*R)/R;
      return Vector2(W*X,W*Y);
   }
   Vector2 acceleration(const Vector2&pos){
      const double X=pos.x;
      const double Y=pos.y;
      const double R=sqrt(X*X+Y*Y);
      base();
      const double W=-pi*pi*c*c*A*C*(c3*R*R*R+c2*R*R+c1*R)/R;
      return Vector2(W*X,W*Y);
   }
   Matrix2 defGrad(const Vector2&pos){
      const double X=pos.x;
      const double Y=pos.y;
      const double R=sqrt(X*X+Y*Y);
      const double cos2=X*X/R/R;
      base();
      const double diag=A*C*X*Y/R*(2.*c3*R+c2);
      return Matrix2(1.+A*C*( 2.*cos2*c3*R*R+cos2*c2*R+   c3*R*R+   c2*R+c1),diag,diag,
                     1.+A*C*(-2.*cos2*c3*R*R-cos2*c2*R+3.*c3*R*R+2.*c2*R+c1));
   }
   Vector2 bodyAccel(const Vector2&pos){
      const double X=pos.x;
      const double Y=pos.y;
      const double R=sqrt(X*X+Y*Y);
      base();

      const double T=A*C;
      const double t1 = pi*pi;
      const double t3 = 1./pch.dens;
      const double t4 = t1*E*t3;
      const double t5 = R*R;
      const double t6 = t5*R;
      const double t8 = c2*t5;
      const double t9 = c1*R;
      const double t11 = T*(c3*t6+t8+t9);
      const double t12 = X/R;
      const double t18 = c2*c2;
      const double t22 = T*c2;
      const double t25 = c3*R;
      const double t27 = T*T;
      const double t28 = t27*T;
      const double t29 = t28*c3;
      const double t34 = c3*c3;
      const double t35 = t28*t34;
      const double t36 = t5*t5;
      const double t41 = c1*c1;
      const double t46 = t27*c3;
      const double t50 = t18*c2;
      const double t54 = t34*t34;
      const double t59 = t27*c2;
      const double t62 = t34*c3;
      const double t63 = t36*R;
      const double t67 = t18*t18;
      const double t72 = t41*c1;
      const double t78 = 6.0*c2+19.0*T*t18*R+12.0*t22*c1+16.0*t25+130.0*t29*t6*c1*t18+221.0*
t35*t36*c1*c2+61.0*t29*t5*t41*c2+122.0*t46*t8*c1+24.0*t27*t50*t5+72.0*t28*t54*
t36*t6+9.0*t59*t41+120.0*t62*t63*t27+12.0*t28*t67*t6+3.0*t28*c2*t72+68.0*T*t34*
t6;
      const double t79 = t27*t41;
      const double t83 = t27*c1;
      const double t86 = T*c3;
      const double t89 = t27*t18;
      const double t94 = t27*t34;
      const double t101 = t28*t62;
      const double t129 = 24.0*t25*t79+112.0*t34*t6*t83+32.0*t86*t9+30.0*t89*t9+76.0*t86*t8+
221.0*t94*t36*c2+191.0*t35*t63*t18+120.0*t101*t63*c1+8.0*t29*R*t72+56.0*t35*t6*
t41+130.0*t46*t6*t18+15.0*t28*t18*t41*R+80.0*t29*t36*t50+195.0*t101*t36*t5*c2+
24.0*t28*t50*c1*t5;
      const double t130 = t78+t129;
      const double t132 = T*c1;
      const double t133 = t22*R;
      const double t135 = t86*t5;
      const double t156 = 1./(t132+1.0+2.0*t133+3.0*t135)/(1.0+2.0*t132+3.0*t59*t9+5.0*t46*t6
*c2+3.0*t133+2.0*t89*t5+4.0*t83*c3*t5+4.0*t135+t79+3.0*t94*t36);
      const double t161 = Y/R;

      return Vector2(-t4*t11*t12-t3*t12*T*E*t130*t156/2.,
                     -t4*t11*t161-t3*t130*E*T*t161*t156/2.);
   }
};
static exact*ex;

void history(patch&pch,ostream&os,int&fc){ // colored splot
   if(single){
      fc+=1;
      os<<'#'<<fc<<'\n';
      for(int i=0;i<pch.Npart();i+=1){
         //os<<pch.px[i]<<'\t'<<radius(pch.px[i]-pch.pX[i],ex->displace(pch.pX[i]))<<'\t'<<pch.param1/double(pch.I)<<'\n';
         //os<<pch.px[i]<<'\t'<<magnitude(pch.pv[i])<<'\t'<<pch.param1/double(pch.I)<<'\n';
         //Matrix2 s;double J;getStress(pch,s,J,i);os<<pch.px[i]<<'\t'<<vonMises(s)<<'\t'<<pch.param1/double(pch.I)<<'\n';

         Matrix2 s;zeroPois(pch,pch.pF[i],s);
         const Vector2 u=unit(pch.px[i]);
         const Matrix2 R (u.x, u.y,-u.y,u.x);
         s=R.inner(s).inner(trans(R));
         //os<<pch.pX[i]+ex->displace(pch.pX[i])<<'\t'<<s.yy<<'\t'<<pch.param1/double(pch.I)<<'\n';
         os<<pch.px[i]<<'\t'<<(pch.pX[i].y<pch.pX[i].x?s.yy:2.*s.xx)<<'\t'<<pch.partSize/double(pch.I)<<'\n';

         //os<<pch.px[i]<<'\t'<<magnitude(pch.pvI[i])<<'\t'<<pch.param1/double(pch.I)<<'\n';
      }
      /* for(int i=0;i<pch.Nnode();++i){
         //if(pch.inRegion(pch.gx[i]))
         os<<pch.gx[i]<<'\t'<<magnitude(pch.ga[i])<<'\t'<<pch.param1/double(pch.I)<<'\n';
      } */
      os<<"\n\n";
   }
   globalLiError=0.;
   globalL2Error=0.;
   globalL1Error=0.;
   for(int i=0;i<pch.Npart();i+=1){
      const double err=radius(pch.px[i]-pch.pX[i],ex->displace(pch.pX[i]));
      //const double err=radius(pch.pxI[i],ex->velocity(pch.pX[i]));
      globalL1Error+=abs(err);
      globalL2Error+=err*err;
      globalLiError=(err>globalLiError?err:globalLiError);
   }
   globalL1Error=globalL1Error/pch.Npart();
   globalL2Error=sqrt(globalL2Error/pch.Npart());
}

void timeIntSC::applyGridBC(patch&pch){
   for(int i=0;i<pch.Nnode();i+=1){
      if(pch.gx[i].x<machTol){
         pch.gv[i].x=0.;
         pch.ga[i].x=0.;
      }
      if(pch.gx[i].y<machTol){
         pch.gv[i].y=0.;
         pch.ga[i].y=0.;
      }
   }
   for(int i=0;i<pch.Nnode();i+=1){
      if(pch.gx[i].x<-machTol){
         pch.gv[i].x=-pch.gv[i+2].x;
         pch.ga[i].x=-pch.ga[i+2].x;
         pch.gv[i].y= pch.gv[i+2].y;
         pch.ga[i].y= pch.ga[i+2].y;
      }
      if(pch.gx[i].y<-machTol){
         pch.gv[i].x= pch.gv[i+2*pch.I].x;
         pch.ga[i].x= pch.ga[i+2*pch.I].x;
         pch.gv[i].y=-pch.gv[i+2*pch.I].y;
         pch.ga[i].y=-pch.ga[i+2*pch.I].y;
      }
   }
}



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
      patch pch(2*Nc,2*Nc,0.,0.,2.,2.,shp.Nghost(),1.);
      timeIntSC&ti=makeTimeInt(tmistr,shp,pch);
      pch.ppe=2.;                      am["ppe"  ]>>pch.ppe;
      pch.load=.1;                     am["load" ]>>pch.load;
      pch.Ymod=1e4;
      pch.dens=1.;
      pch.vwav=sqrt(pch.Ymod/pch.dens);
      pch.pois=.0;
      pch.damp=.0;
      pch.partSize=32.;                am["p1"   ]>>pch.partSize;
      pch.dt=min(pch.dx,pch.dy)*CFL/sqrt(pch.Ymod/pch.dens);
      //fillAnnulusJCP(pch,Vector2(0.,0.),.5,1.,pch.ppe,shp.Nsupport());
      fillAnnulusRegular(pch,Vector2(0.,0.),.5,1.,pch.ppe,shp.Nsupport());
      //fillAnnulusRadial(pch,Vector2(0.,0.),.5,1.,shp.Nsupport());
      //fillRefinedSurface(pch,inRing,int(pch.ppe),2*int(pch.ppe),shp.Nsupport());
         double vsum=0.;
         for(int i=0;i<pch.Npart();++i)vsum+=pch.pV[i];
         const double asum=vsum/pch.thick;
         cerr.precision(14);
         cerr<<"area ratio: "<<asum/(.1875*pi)<<tab;
         ex=new exact(pch); // do after all patch params assigned
      pch.elaps=.0/pch.vwav;
      for(int i=0;i<pch.Npart();i+=1){
         pch.px[i]=pch.pX[i]+ex->displace(pch.pX[i]);
         pch.pv[i]=ex->velocity(pch.pX[i]);
         pch.pF[i]=ex->defGrad(pch.pX[i]);
      }
      shp.updateContribList(pch);
// main loop
      ofstream hf("history.xls");
      int frameCount=0;
      history(pch,hf,frameCount);
      double Li=0.,L2=0.,L1=0.;
      try{
         while(pch.elaps<2./pch.vwav){
         //while(pch.incCount<10){
            for(int i=0;i<pch.Npart();++i)pch.pfe[i]=(tmistr=="MLS"?1.:pch.pm[i])*ex->bodyAccel(pch.pX[i]);
            ti.advance(pch);
            pch.elaps+=pch.dt;
            pch.incCount+=1;
            history(pch,hf,frameCount);
            L1=(globalL1Error>L1?globalL1Error:L1);
            L2=(globalL2Error>L2?globalL2Error:L2);
            Li=(globalLiError>Li?globalLiError:Li);
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
         cerr<<"Wall time:"<<double(clock()-t)/CPS<<tab<<"L-inf:"<<Li<<endl;
         cout<<frameCount;
      }
      else{
         cerr<<"Wall time:"<<double(clock()-t)/CPS<<tab<<"L-inf:"<<Li<<endl;
         cout<<pch.Npart()<<'\t'
             <<pch.I*pch.J<<'\t'
             <<pch.incCount<<'\t'
             <<double(clock()-t)/CPS<<'\t'
             <<L1<<'\t'
             <<L2<<'\t'
             <<Li;
      }
   }catch(...){cerr<<"Non-loop Exception"<<endl;}
}


