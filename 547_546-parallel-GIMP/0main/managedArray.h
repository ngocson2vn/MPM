// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef MANAGED_ARRAY_H
#define MANAGED_ARRAY_H

#include<tensor3.h>
using namespace std;

inline MPI_Datatype*MPItype(realT){return &MPIrealT;} // because MPI_DOUBLE is not an lvalue
inline MPI_Datatype*MPItype(vec3){return &MPIvec3;}
inline MPI_Datatype*MPItype(mat3){return &MPImat3;}
typedef void* voidPtr;

// operations required of all arrays in an arrayGroup
struct arraySC{
   virtual void resize(size_t)=0;
   virtual indexT size()=0;
   virtual void toBuff(indexT i,voidPtr&)=0;
   virtual void fromBuff(voidPtr&)=0;
   virtual ~arraySC(){}
};

// each array of type T uses the operation from its underlying vector
template<typename T>struct managedArray:public arraySC,public vector<T>{
   void resize(size_t s){vector<T>::resize(s);}
   indexT size(){return indexT(vector<T>::size());}
   void toBuff(indexT p,voidPtr&vd){
      T*ptr=static_cast<T*>(vd);     // cast to element type
      *ptr=vector<T>::operator[](p); // copy element to buffer
      ++ptr;                         // accounts for element size
      vd=ptr;                        // vd now incremented by size of element
      vector<T>::operator[](p)=vector<T>::back(); // copy over element with last
      vector<T>::pop_back();         // delete last
   }
   void fromBuff(void*&vd){
      T*ptr=static_cast<T*>(vd);     // cast to element type
      vector<T>::push_back(*ptr);    // add element to end of array
      ++ptr;                         // accounts for element size
      vd=ptr;                        // vd now incremented by size of element
   }
};

// a heterogeneous group of arrays, accessed via virtual functions
class arrayGroup{
   indexT eachSize;
   vector<MPI_Datatype*>mpiTypes; // only want one of each type, so we point to it, instead of copying
   vector<arraySC*>ag;
public:
   arrayGroup(){eachSize=0;}
   void resizeAll(const indexT s){
      eachSize=s;
      for(indexT i=0;i<ag.size();++i)ag[i]->resize(s);
   }
   void verifySize(){
      for(indexT i=0;i<ag.size();++i){
         assert(ag[i]->size() == eachSize);
      }
   }
   indexT size(){return eachSize;}
   void addArray(arraySC*p,MPI_Datatype*m){
      mpiTypes.push_back(m);
      ag.push_back(p);
      p->resize(eachSize);
   }
   void makeMPIbundleType(MPI_Datatype*bundleType){
      //MPI_Type_dup(MPI_DOUBLE,&MPIrealT);
      //MPI_Type_dup(MPI_UNSIGNED,&MPIindexT);
      //MPI_Type_contiguous(3,MPI_DOUBLE,&MPIvec3);
      //MPI_Type_commit(&MPIvec3);
      //MPI_Type_contiguous(9,MPI_DOUBLE,&MPImat3);
      //MPI_Type_commit(&MPImat3);
      int count=static_cast<int>(mpiTypes.size());
      int*blocklen=new int[count];
      for(int c=0;c<count;++c)blocklen[c]=1;
      MPI_Aint*disp=new MPI_Aint[count];
      disp[0]=0;
      for(int c=1;c<count;++c){
         int typeSize;
         MPI_Type_size(*mpiTypes[c-1],&typeSize);
         disp[c]=disp[c-1]+typeSize;
      }
      MPI_Datatype*types=new MPI_Datatype[count];
      for(int c=0;c<count;++c)types[c]=*mpiTypes[c]; // now we copy the types
      MPI_Type_create_struct(count,blocklen,disp,types,bundleType);
      MPI_Type_commit(bundleType);
   }
   void partToBuff(indexT p,voidPtr&vd){
      for(indexT i=0;i<ag.size();++i)ag[i]->toBuff(p,vd);
      --eachSize;
   }
   void partFromBuff(voidPtr&vd){
      for(indexT i=0;i<ag.size();++i)ag[i]->fromBuff(vd);
      ++eachSize;
   }
};

// two global array groups; arrays may be added anytime, from anywhere
extern arrayGroup globalPartArrays;
extern arrayGroup globalNodeArrays;

// convenience classes for defining arrays that belong to a managed group
template<typename T>
struct partArray:public managedArray<T>{
   partArray(){
      T dummy;
      dummy=0.;
      globalPartArrays.addArray(this,MPItype(dummy));
   }
};

template<typename T>
struct nodeArray:public managedArray<T>{
   nodeArray(){
      T dummy;
      dummy=0.;
      globalNodeArrays.addArray(this,MPItype(dummy));
   }
};

// convenience functions for array group sizes
inline indexT Npart(){return globalPartArrays.size();}
inline indexT Nnode(){return globalNodeArrays.size();}

// unused diagnostic
inline indexT Nmask(const nodeArray<indexT>&mask){
   indexT count=0;
   for(indexT n=0;n<Nnode();++n)if(mask[n]!=0)++count;
   return count;
}


#endif
