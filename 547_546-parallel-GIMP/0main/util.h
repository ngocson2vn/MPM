// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef UTIL_H
#define UTIL_H

#include<exception>
#include<limits>
#include<iomanip>
#include<fstream>
#include<vector>
#include<map>
#include<sstream>
#include<complex>
#include<mpi.h> // must come before BUILTIN_DOUBLE_PROHIBITED
using namespace std;

typedef double realT;
typedef unsigned indexT;
#define float BUILTIN_FLOAT_PROHIBITED
#define double BUILTIN_DOUBLE_PROHIBITED
#define unsigned BUILTIN_UNSIGNED_PROHIBITED
#define short BUILTIN_SHORT_PROHIBITED

struct stop:public vector<string>{
   stop(const char*a){push_back(a);}
   void operator+=(const char*a){push_back(a);}
};

const realT pi=3.1415926535897932384626433832795028841971694;
const realT machEps=numeric_limits<realT>::epsilon();
const realT oneEps=1.-machEps;
const realT tiny=numeric_limits<realT>::min();
const realT huge=numeric_limits<realT>::max();

template<typename T>inline T sgn(T a){return(a<0.?-1.:1.);}
template<bool>struct Assert;
template<>struct Assert<true>{}; // Assert<(1==2)>(); <-- will not compile

extern MPI_Comm cartComm;
extern int cartRank,cartSize;
extern MPI_Datatype MPIindexT,MPIrealT,MPIvec3,MPImat3,MPIparticle;
extern MPI_Op MPI_DDadd;
void here(const int h);
void crash(stop&sp);
void commCheck();
void waitCheckMPI(MPI_Request*,MPI_Status*);

inline void DDadd(void*invec,void*inoutvec,int*,MPI_Datatype*){
   // Uses the double-double precision of
   // Yun He and Chris Ding, 12/06/99, NERSC
   complex<realT>*X=static_cast<complex<realT>*>(invec);
   complex<realT>*S=static_cast<complex<realT>*>(inoutvec);
   const realT t1=X->real()+S->real();
   const realT e=t1-X->real();
   const realT t2=((S->real()-e)+(X->real()-(t1-e)))+X->imag()+S->imag();
   *S=complex<realT>(t1+t2,t2-((t1+t2)-t1));
}

// custom stream operators for tab and non-flushed newline
template<class charT,class traits>basic_ostream<charT,traits>&tab(basic_ostream<charT,traits>&os){os<<'\t'<<setw(14);return os;}
template<class charT,class traits>basic_ostream<charT,traits>&nwl(basic_ostream<charT,traits>&os){os<<'\n';return os;}

template<typename T>string toStr(const T&val){ // universal conversion to string
   ostringstream osval;
   osval.precision(14);
   osval<<val;
   return osval.str();
}

struct timeTag{const realT sec;timeTag(const realT tm):sec(tm){}};
inline std::ostream& operator<<(std::ostream&os,const timeTag tr){
   const int w=4;
   const realT&t=tr.sec;
   if(t>3600.){os<<setw(w)<<t/3600.<<" hr  ";return os;}
   if(t>60.  ){os<<setw(w)<<t/60.  <<" min ";return os;}
               os<<setw(w)<<t      <<" sec ";return os;
}

//////// command-line args /////////////

struct argRec{
   const string val;
   mutable indexT Naccess;
   argRec(const string&a):val(a){Naccess=0;}
};

class CLargMap:private map<const string,const argRec>{
   template<typename T>friend T comLineArg(const string&tag,T val);
public:
   void init(int argc,char**argv){
      for(int arg=1;arg<argc;++arg){
         string argstr(argv[arg]);
         if(argstr[0]=='-'){ // skipping posix leading dash args
            argstr=argv[arg+1];
            if(argstr.find('=')==string::npos)arg+=1; // possibly skip next arg too
            continue;
         }
         for(indexT j=0;j<argstr.size();++j)if(argstr[j]=='~')argstr[j]=' ';
         if(argstr.find('=')==string::npos)throw stop("invalid stdin - did you forget to put '~' between vector components?");
         const string tag(argstr.substr(0,argstr.find('=')));    // exception if = not found
         const string val(argstr.substr(tag.size()+1));
         if(find(tag)==end())insert(make_pair(tag,argRec(val)));
         else throw stop("duplicate command line tag encountered");
      }
   }
   void reportUnusedTags(){
      for(CLargMap::const_iterator it=begin();it!=end();++it){
         if(it->second.Naccess<1)cerr<<"Warning, unused command line argument: "<<it->first<<" = "<<it->second.val<<endl;
      }
   }
};
extern CLargMap globalArgMap;

//////// tag=val reports /////////////

class paramReport:public map<string,string>{
   ofstream reportFile;
   string reportName;
public:
   void init(){
      reportName=comLineArg("testName",string("run"))+string(".report");
      write();
   }
   template<typename T>void param(const string&a,const T&b){operator[](a)=toStr(b);}
   void write(){
      if(cartRank==0){
         reportFile.open(reportName.c_str());
         reportFile.precision(14);
         for(map<string,string>::iterator it=begin();it!=end();++it){reportFile<<it->first<<'='<<it->second<<nwl;}
         reportFile.close();
      }
   }
   void writeFinal(){
      for(int rank=0;rank<cartSize;++rank){
         if(rank==cartRank){
            if(rank==0)reportFile.open(reportName.c_str(),ios::out);
            else       reportFile.open(reportName.c_str(),ios::app);
            reportFile.precision(14);
            reportFile<<"cartRank_"<<cartRank<<"=_____________________________"<<nwl;
            for(map<string,string>::iterator it=begin();it!=end();++it){reportFile<<it->first<<'='<<it->second<<nwl;}
            reportFile.close();
         }
         MPI_Barrier(cartComm);
      }
   }
};

extern paramReport report;

template<typename T>
T comLineArg(const string&tag,T val){
   CLargMap::const_iterator it=globalArgMap.find(tag);
   if(it==globalArgMap.end()){
      report.param(tag,toStr(val));
      return val; // tag not found, return default val
   }
   istringstream sval(it->second.val);
   sval>>val;
   ++it->second.Naccess;      // increment the access counter
   report.param(tag,toStr(val));
   return val;                // tag found, return command line val
}

// Convenience class for a geometric series of load factors
// Useful for non-linear problems such as pull-in
class geomStepSeries{
   realT dt;
   realT beta;
   int N;
   bool alreadyInit;
public:
   geomStepSeries(){alreadyInit=false;}
   void init(realT lastFactor,int totSteps){
      if(alreadyInit)throw stop("geomStepSeries may only be initialized once!");
      dt=lastFactor;
      N=totSteps;
      if(totSteps>1 && abs(1./realT(totSteps)-lastFactor)>machEps){
         beta=1.1;
         realT fcurr=1.-dt*(pow(beta,N)-1.)/(beta-1.);
         realT fprev=0.;
         for(int c=0;c<2 || abs(fcurr)<fprev;++c){
            fprev=abs(fcurr);
            const realT df=dt*((pow(beta,N)-1.)/((beta-1.)*(beta-1.))-(pow(beta,N)*realT(N))/(beta*(beta-1.)));
            beta-=fcurr/df;
            fcurr=1.-dt*(pow(beta,N)-1.)/(beta-1.);
            cerr<<"geomStepSeries: Newton-Raphson: f / beta: "<<fcurr<<" / "<<beta<<endl;
         }
         cerr<<"geomStepSeries: beta: "<<beta<<endl;
         cerr<<"geomStepSeries: first factor: "<<factorI(0)<<endl;
      }
      else beta=1.;
      alreadyInit=true;
   }
   realT factorI(const int I)const{
      if(beta==1.)return realT(I+1)*dt;
      return 1.-dt*(pow(beta,N-(I+1))-1.)/(beta-1.); // big to small
      //return dt*(pow(beta,(I+1))-1.)/(beta-1.);    // small to big (should check before use)
   }
};

#endif
