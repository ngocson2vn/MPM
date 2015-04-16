
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

// This is the constitutive model for this problem.
double planeStrainNeoHookean(const patch&pch,const Matrix2&F,Matrix2&S){
   const double v=pch.pois;
   const double l=pch.Ymod*v/((1.+v)*(1.-2.*v));
   const double m=.5*pch.Ymod/(1.+v);
   const double Ja=det(F);
   S =I2(l*log(Ja)/Ja);
   S+=m/Ja*(F.inner(trans(F))-I2());
   return Ja;
}
// This fulfills the virtual function requirement for the timeIntSC
// super class.  It gets called at least once per particle per time
// step, which is bad.  I do it so the user can specify the 
// constitutive model, but hopefully there is a better way
double timeIntSC::getStress(const patch&pch,const Matrix2&F,Matrix2&S){return planeStrainNeoHookean(pch,F,S);}

// For this problem there is an exact manufactured solution.
// This allows us to measure the convergence rate of the method.
struct exact{
   const patch&pch;
   double X,Y,t,S,C,v,c,A,Cl,Cm;
   exact(const patch&p):pch(p){}
   void base(){
      t=pch.elaps;
      c=pch.vwav;
      S=sin(c*pi*t);
      C=cos(c*pi*t);
      v=pch.pois;
      A=pch.load;
      Cl=v/((1.+v)*(1.-2.*v));
      Cm=.5/(1.+v);
   }
   Vector2 getX(const Vector2&n){
      base();
      double X=n.x;
      while(true){
         const double f=X+A*C*sin(pi*X)-n.x;
         if(abs(f)<machTol)break;
         X-=f/(1.+A*C*pi*cos(pi*X));
      }
      double Y=n.y;
      while(true){
         const double f=Y+A*S*sin(pi*Y)-n.y;
         if(abs(f)<machTol)break;
         Y-=f/(1.+A*S*pi*cos(pi*Y));
      }
      return Vector2(X,Y);
   }
   Vector2 displace(const Vector2&pos){
      X=pos.x;
      Y=pos.y;
      base();
      return A*Vector2(sin(pi*X)*C,sin(pi*Y)*S);
   }
   Vector2 velocity(const Vector2&pos){
      X=pos.x;
      Y=pos.y;
      base();
      return A*c*pi*Vector2(sin(pi*X)*S,sin(pi*Y)*C);
   }
   Vector2 accel(const Vector2&pos){
      X=pos.x;
      Y=pos.y;
      base();
      return -A*c*c*pi*pi*Vector2(sin(pi*X)*C,sin(pi*Y)*S);
   }
   Matrix2 defGrad(const Vector2&pos){
      X=pos.x;
      Y=pos.y;
      base();
      return Matrix2(1.+A*pi*cos(pi*X)*C,0.,0.,1.+A*pi*cos(pi*Y)*S);
   }
   Vector2 bodyAccel(const Vector2&pos){
      X=pos.x;
      Y=pos.y;
      const double Fx=1.+A*pi*cos(pi*X)*C;
      const double Fy=1.+A*pi*cos(pi*Y)*S;
      const Vector2 u=A*Vector2(sin(pi*X)*C,sin(pi*Y)*S);
      base();
      return pch.Ymod/pch.dens*pi*pi*u*Vector2((1.-log(Fx*Fy))*Cl/Fx/Fx+(1.+1./Fx/Fx)*Cm-1.,
                                               (1.-log(Fx*Fy))*Cl/Fy/Fy+(1.+1./Fy/Fy)*Cm-1.);
   }
};
static exact*ex;

// all-around "print stuff" function
void history(patch&pch,ostream&os,int&fc){ // colored splot
   // If its a single run, we write data to disk.  Otherwise we
   // only collect measures of error.
   if(single){
      fc+=1;
      os<<'#'<<fc<<'\n';
      for(int i=0;i<pch.Npart();++i){
         os<<pch.px[i]<<'\t'<<pch.pJ[i]<<'\t'<<pch.partSize/double(pch.I)<<'\n';
      }
      os<<"\n\n";
   }
   // Here we measure error as the magnitude of the vector of
   // displacement error: ||u_actual - u_exact||
   globalLiError=0.;
   globalL2Error=0.;
   globalL1Error=0.;
   for(int i=0;i<pch.Npart();++i){
      const double err=radius(pch.px[i]-pch.pX[i],ex->displace(pch.pX[i]));
      globalL1Error+=abs(err);
      globalL2Error+=err*err;
      globalLiError=(err>globalLiError?err:globalLiError);
   }
   globalL1Error=globalL1Error/pch.Npart();
   globalL2Error=sqrt(globalL2Error/pch.Npart());
}

// Another virtual function.  This one is not so bad,
// because it is only called once for the whole loop
void timeIntSC::applyGridBC(patch&pch){
   for(int i=0;i<pch.Nnode();i+=1){
      if(pch.gx[i].x<machTol||pch.gx[i].x>1.-machTol){
         pch.gv[i].x=0.;
         pch.ga[i].x=0.;
      }
      if(pch.gx[i].y<machTol||pch.gx[i].y>1.-machTol){
         pch.gv[i].y=0.;
         pch.ga[i].y=0.;
      }
   }
}

int main(int argc,char**argv){
   // Enclose the entire program in a try block so we have a chance
   // of getting information in the event of a crash
   try{

      // typical calls to measure the wall time
      const double CPS=double(CLOCKS_PER_SEC);
      clock_t t=clock();
      
      // assemble command line parameter map
      argMap am(argc,argv); // assemble command line parameter map
      
      // Guesses are made for all parameters and the guesses are
      // overwritten by the command line argument, if the appropriate
      // tag is found.  Most of the tags below are described in the README.
      // The format is:
      // guess a value here            if the tag is found, then insert the command line value into the variable
      int Nc=20;                       am["Ncell"]>>Nc;
      string shpstr="GIMP";            am["shape"]>>shpstr;
      string tmistr="cen";             am["tInt" ]>>tmistr;
      double CFL=.2;                   am["CFL"  ]>>CFL;
                                       am["single"]>>single;
      
      // Based on the shape string a new instance of shape functions and operators
      // is created and a reference to it is returned
      shapeSC&shp=makeShape(shpstr);
      
      // Omnibus patch constructor.  Arguments are:
      // 1. number of cells in x-dir
      // 2. number of cells in y-dir
      // 3. left boundary of domain
      // 4. lower boundary of domain
      // 5. right boundary of domain
      // 6. upper boundary of domain
      // 7. number of layers of ghost cells
      // 8. assumed thickness of material
      patch pch(Nc,Nc,0.,0.,1.,1.,shp.Nghost(),1.);
      
      // A new instance of the desired time-stepping class is created and
      // a reference to it is returned.
      timeIntSC&ti=makeTimeInt(tmistr,shp,pch);
      
      // make a new exact solution class, if one is defined
      ex=new exact(pch);
      
      // initialize various physical parameters
      pch.ppe=2.;                      am["ppe"  ]>>pch.ppe;
      pch.load=.1;                     am["load" ]>>pch.load;
      pch.Ymod=1e4;
      pch.dens=1.;
      pch.vwav=sqrt(pch.Ymod/pch.dens);
      pch.pois=.3;
      pch.damp=.0;
      pch.partSize=32.;                am["p1"   ]>>pch.partSize; // gnuplot point size
      
      // This is how we calculate the time step size based on the
      // characteristic cell size, the material wave speed, and the desired CFL
      pch.dt=min(pch.dx,pch.dy)*CFL/sqrt(pch.Ymod/pch.dens);
      
      // This function fills a rectangle with particles.  We may call this
      // and the other functions in io.cpp to form various shapes.
      // It is impossible to predict the number of particles that will be
      // needed in the simulation until after the shapes are formed.
      // If shapes overlap then extraneous particles must be trimmed.
      fillRectangleRegular(pch,Vector2(0.,0.),Vector2(1.,1.),pch.ppe,shp.Nsupport());
      
      // We initialize the simulation time to zero.  But for the cyclical
      // manufactured solutions we might want to start at pi/2 or whatever
      pch.elaps=.0/pch.vwav;
      
      // Initializing particles to the exact values of the manufactured solution
      // if there is one, or to some known initial state
      for(int i=0;i<pch.Npart();i+=1){
         pch.px[i]=pch.pX[i]+ex->displace(pch.pX[i]);
         pch.pv[i]=ex->velocity(pch.pX[i]);
         pch.pF[i]=ex->defGrad(pch.pX[i]);
      }
      
      // Figure out what cell surrounds each particle and pre-calculate the
      // shapes and gradients weights for each particle to each node
      shp.updateContribList(pch);

      ofstream hf("history.xls"); // an output file
      int frameCount=0;           // not every time step outputs data
      double Li=0.,L2=0.,L1=0.;   // measures of error
      
      // Another try block with more detailed handling of crashes
      try{

         // 2/vwav would be the wave transit time for a rod of length 2
         // this is just a convenient characteristic time scale for these problems
         while(pch.elaps<2./pch.vwav){
         //while(pch.incCount<1){

            // impose known body forces on particles
            for(int i=0;i<pch.Npart();++i)pch.pfe[i]=pch.pm[i]*ex->bodyAccel(pch.pX[i]);
            
            // advance forward one time step (this is the core of the algorithm)
            ti.advance(pch);

            pch.elaps+=pch.dt;          // advanced the elapsed time
            pch.incCount+=1;            // advance the iteration count                          
            history(pch,hf,frameCount); // write data to disk, and get errors for this time step

            // continually look for the highest error
            L1=(globalL1Error>L1?globalL1Error:L1);  // L-1 norm
            L2=(globalL2Error>L2?globalL2Error:L2);  // L-2 norm
            Li=(globalLiError>Li?globalLiError:Li);  // L-infinity norm
            
            // updating the contribution nodes and weights
            // I do this last instead of first because this is where the code
            // is designed to crash if particles have strayed out of bounds
            shp.updateContribList(pch);

            cerr<<'~'; // show the user we made a step
         }
      }
      
      // "Defined Exceptions" are purposely thrown by the code
      catch(const char*s){cerr<<"Defined Exception:"<<s<<endl;L1=L2=Li=1.;}
      
      // "Standard Exceptions" are thrown by the C++ standard library
      catch(const std::exception&error){cerr<<"Standard Exception:"<<error.what()<<endl;}
      
      // "Unknown Exceptions" are just that
      catch(...){cerr<<"Unknown Exception"<<endl;}
      cerr<<'\n';

      // A "single" run solves the problem at one mesh size and sends the number
      // of time steps written to disk (not necessarily all of them) to stdout
      if(single){
         cerr<<"Wall time:"<<double(clock()-t)/CPS<<tab<<"L-inf:"<<Li<<endl;
         cout<<frameCount;
      }
      
      // Non-single means we're running the program over a series of mesh sizes
      // in order to make a convergence plot.  So we don't want to write data
      // to disk in the history() function and we do want to send various
      // performance measures to stdout, where python can receive them.
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

   // end of the outer-most try block
   }catch(...){cerr<<"Non-loop Exception"<<endl;}
}


