// Copyright (c) 2009,2011, Philip Wallstedt.  All rights reserved.
// This software is covered by the license and liability disclaimer in "license.txt"

#ifndef TENSOR_H
#define TENSOR_H

#include<iostream>  // printing vec and mat
#include<iomanip>
#include<cassert>   // unit
#include<cmath>     // sqrt
#include"util.h"

using namespace std;

struct mat3;
struct vec3{
   realT x,y,z;
   vec3(){x=0.;y=0.;z=0.;}
   vec3(const realT a,const realT b,const realT c){x=a;y=b;z=c;}
   vec3(const vec3&v){x=v.x;y=v.y;z=v.z;}
   realT&operator[](const indexT i){
      switch(i){
         case 0:return x;
         case 1:return y;
         case 2:return z;
         default:assert(false);
      }
   }
   const realT&operator[](const indexT i)const{
      switch(i){
         case 0:return x;
         case 1:return y;
         case 2:return z;
         default:assert(false);
      }
   }
   void roll(){const realT x0=x;x=y;y=z;z=x0;}
   vec3&operator=(const realT&a){x=a;y=a;z=a;return*this;}
   vec3&operator=(const vec3&v){x=v.x;y=v.y;z=v.z;return*this;}
   vec3&operator+=(const vec3&v){x+=v.x;y+=v.y;z+=v.z;return*this;}
   vec3&operator-=(const vec3&v){x-=v.x;y-=v.y;z-=v.z;return*this;}
   vec3&operator*=(const vec3&v){x*=v.x;y*=v.y;z*=v.z;return*this;}
   vec3&operator/=(const vec3&v){x/=v.x;y/=v.y;z/=v.z;return*this;}
   vec3&operator+=(const realT&v){x+=v;y+=v;z+=v;return*this;}
   vec3&operator-=(const realT&v){x-=v;y-=v;z-=v;return*this;}
   vec3&operator*=(const realT&v){x*=v;y*=v;z*=v;return*this;}
   vec3&operator/=(const realT&v){x/=v;y/=v;z/=v;return*this;}
   const vec3 operator-(            )const{return vec3(-x,-y,-z);}
   const vec3 operator+(const vec3&v)const{return vec3(x+v.x,y+v.y,z+v.z);}
   const vec3 operator-(const vec3&v)const{return vec3(x-v.x,y-v.y,z-v.z);}
   const vec3 operator*(const vec3&v)const{return vec3(x*v.x,y*v.y,z*v.z);}
   const vec3 operator/(const vec3&v)const{return vec3(x/v.x,y/v.y,z/v.z);}
   const vec3 operator+(const realT&a)const{return vec3(x+a,y+a,z+a);}
   const vec3 operator-(const realT&a)const{return vec3(x-a,y-a,z-a);}
   const vec3 operator*(const realT&a)const{return vec3(x*a,y*a,z*a);}
   const vec3 operator/(const realT&a)const{return vec3(x/a,y/a,z/a);}
   bool operator==(const vec3&v)const{return x==v.x&&y==v.y&&z==v.z;}
   bool operator==(const realT &v)const{return x==v&&y==v&&z==v;}
   bool operator!=(const vec3&v)const{return x!=v.x||y!=v.y||z!=v.z;}
   bool operator!=(const realT &v)const{return x!=v||y!=v||z!=v;}
   bool operator<=(const vec3&v)const{return x<=v.x&&y<=v.y&&z<=v.z;}
   bool operator< (const vec3&v)const{return x< v.x&&y< v.y&&z< v.z;}
   inline realT  inner(const vec3&v)const;
   inline vec3 inner(const mat3&m)const;
   inline mat3 outer(const vec3&v)const;
};
inline const vec3 operator+(const realT&a,const vec3&v){return vec3(a+v.x,a+v.y,a+v.z);}
inline const vec3 operator-(const realT&a,const vec3&v){return vec3(a-v.x,a-v.y,a-v.z);}
inline const vec3 operator*(const realT&a,const vec3&v){return vec3(a*v.x,a*v.y,a*v.z);}
inline const vec3 operator/(const realT&a,const vec3&v){return vec3(a/v.x,a/v.y,a/v.z);}

struct mat3{
   realT xx,xy,xz,yx,yy,yz,zx,zy,zz;
   mat3(){xx=0.;xy=0.;xz=0.;yx=0.;yy=0.;yz=0.;zx=0.;zy=0.;zz=0.;}
   mat3(realT a){xx=a;xy=a;xz=a;yx=a;yy=a;yz=a;zx=a;zy=a;zz=a;}
   mat3(const realT a,const realT b,const realT c,
        const realT d,const realT e,const realT f,
        const realT g,const realT h,const realT i){xx=a;xy=b;xz=c;yx=d;yy=e;yz=f;zx=g;zy=h;zz=i;}
   mat3(const mat3&m){xx=m.xx;xy=m.xy;xz=m.xz;yx=m.yx;yy=m.yy;yz=m.yz;zx=m.zx;zy=m.zy;zz=m.zz;}
   realT operator()(const indexT i,const indexT j)const{
      const indexT k=3*i+j;
      switch(k){
         case 0:return xx;
         case 1:return xy;
         case 2:return xz;
         case 3:return yx;
         case 4:return yy;
         case 5:return yz;
         case 6:return zx;
         case 7:return zy;
         case 8:return zz;
         default:assert(false);
      }
   }
   realT&operator()(const indexT i,const indexT j){
      const indexT k=3*i+j;
      switch(k){
         case 0:return xx;
         case 1:return xy;
         case 2:return xz;
         case 3:return yx;
         case 4:return yy;
         case 5:return yz;
         case 6:return zx;
         case 7:return zy;
         case 8:return zz;
         default:assert(false);
      }
   }
   mat3&operator=(const realT&a){xx=a;xy=a;xz=a;yx=a;yy=a;yz=a;zx=a;zy=a;zz=a;return*this;}
   mat3&operator=(const mat3&m){xx=m.xx;xy=m.xy;xz=m.xz;yx=m.yx;yy=m.yy;yz=m.yz;zx=m.zx;zy=m.zy;zz=m.zz;return*this;}
   mat3&operator+=(const mat3&m){xx+=m.xx;xy+=m.xy;xz+=m.xz;yx+=m.yx;yy+=m.yy;yz+=m.yz;zx+=m.zx;zy+=m.zy;zz+=m.zz;return*this;}
   mat3&operator-=(const mat3&m){xx-=m.xx;xy-=m.xy;xz-=m.xz;yx-=m.yx;yy-=m.yy;yz-=m.yz;zx-=m.zx;zy-=m.zy;zz-=m.zz;return*this;}
   mat3&operator*=(const mat3&m){xx*=m.xx;xy*=m.xy;xz*=m.xz;yx*=m.yx;yy*=m.yy;yz*=m.yz;zx*=m.zx;zy*=m.zy;zz*=m.zz;return*this;}
   mat3&operator/=(const mat3&m){xx/=m.xx;xy/=m.xy;xz/=m.xz;yx/=m.yx;yy/=m.yy;yz/=m.yz;zx/=m.zx;zy/=m.zy;zz/=m.zz;return*this;}
   mat3&operator+=(const realT&m){xx+=m;xy+=m;xz+=m;yx+=m;yy+=m;yz+=m;zx+=m;zy+=m;zz+=m;return*this;}
   mat3&operator-=(const realT&m){xx-=m;xy-=m;xz-=m;yx-=m;yy-=m;yz-=m;zx-=m;zy-=m;zz-=m;return*this;}
   mat3&operator*=(const realT&m){xx*=m;xy*=m;xz*=m;yx*=m;yy*=m;yz*=m;zx*=m;zy*=m;zz*=m;return*this;}
   mat3&operator/=(const realT&m){xx/=m;xy/=m;xz/=m;yx/=m;yy/=m;yz/=m;zx/=m;zy/=m;zz/=m;return*this;}
   const mat3 operator-(            )const{return mat3(-xx,-xy,-xz,-yx,-yy,-yz,-zx,-zy,-zz);}
   const mat3 operator+(const mat3&m)const{return mat3(xx+m.xx,xy+m.xy,xz+m.xz,yx+m.yx,yy+m.yy,yz+m.yz,zx+m.zx,zy+m.zy,zz+m.zz);}
   const mat3 operator-(const mat3&m)const{return mat3(xx-m.xx,xy-m.xy,xz-m.xz,yx-m.yx,yy-m.yy,yz-m.yz,zx-m.zx,zy-m.zy,zz-m.zz);}
   const mat3 operator*(const mat3&m)const{return mat3(xx*m.xx,xy*m.xy,xz*m.xz,yx*m.yx,yy*m.yy,yz*m.yz,zx*m.zx,zy*m.zy,zz*m.zz);}
   const mat3 operator/(const mat3&m)const{return mat3(xx/m.xx,xy/m.xy,xz/m.xz,yx/m.yx,yy/m.yy,yz/m.yz,zx/m.zx,zy/m.zy,zz/m.zz);}
   const mat3 operator+(const realT&m)const{return mat3(xx+m,xy+m,xz+m,yx+m,yy+m,yz+m,zx+m,zy+m,zz+m);}
   const mat3 operator-(const realT&m)const{return mat3(xx-m,xy-m,xz-m,yx-m,yy-m,yz-m,zx-m,zy-m,zz-m);}
   const mat3 operator*(const realT&m)const{return mat3(xx*m,xy*m,xz*m,yx*m,yy*m,yz*m,zx*m,zy*m,zz*m);}
   const mat3 operator/(const realT&m)const{return mat3(xx/m,xy/m,xz/m,yx/m,yy/m,yz/m,zx/m,zy/m,zz/m);}
   bool operator==(const mat3&m){return xx==m.xx&&xy==m.xy&&xz==m.xz&&yx==m.yx&&yy==m.yy&&yz==m.yz&&zx==m.zx&&zy==m.zy&&zz==m.zz;}
   bool operator!=(const mat3&m){return xx!=m.xx||xy!=m.xy||xz!=m.xz||yx!=m.yx||yy!=m.yy||yz!=m.yz||zx!=m.zx||zy!=m.zy||zz!=m.zz;}
   inline vec3 inner(const vec3&v)const;
   inline mat3 inner(const mat3&n)const;
};
inline const mat3 operator+(const realT&a,const mat3&m){return mat3(a+m.xx,a+m.xy,a+m.xz,a+m.yx,a+m.yy,a+m.yz,a+m.zx,a+m.zy,a+m.zz);}
inline const mat3 operator-(const realT&a,const mat3&m){return mat3(a-m.xx,a-m.xy,a-m.xz,a-m.yx,a-m.yy,a-m.yz,a-m.zx,a-m.zy,a-m.zz);}
inline const mat3 operator*(const realT&a,const mat3&m){return mat3(a*m.xx,a*m.xy,a*m.xz,a*m.yx,a*m.yy,a*m.yz,a*m.zx,a*m.zy,a*m.zz);}
inline const mat3 operator/(const realT&a,const mat3&m){return mat3(a/m.xx,a/m.xy,a/m.xz,a/m.yx,a/m.yy,a/m.yz,a/m.zx,a/m.zy,a/m.zz);}

//       unary operators                                   //
inline realT len2(const vec3&u){return u.x*u.x+u.y*u.y+u.z*u.z;}
inline realT trace(const mat3&m){return m.xx+m.yy+m.zz;}
inline vec3 diag(const mat3&m){return vec3(m.xx,m.yy,m.zz);}
inline vec3 round(const vec3&v){return vec3(round(v.x),round(v.y),round(v.z));}
inline realT dotDot(const mat3&m,const mat3&n){return m.xx*n.xx+m.xy*n.xy+m.xz*n.xz
                                                     +m.yx*n.yx+m.yy*n.yy+m.yz*n.yz
                                                     +m.zx*n.zx+m.zy*n.zy+m.zz*n.zz;}
inline realT Frob(const mat3&m){return sqrt(dotDot(m,m));}
inline realT det(const mat3&m){return m.xx*m.yy*m.zz+m.xy*m.yz*m.zx+m.xz*m.yx*m.zy-m.zx*m.yy*m.xz-m.zy*m.yz*m.xx-m.zz*m.yx*m.xy;}
inline mat3 I3(const realT f=1.){return mat3(f,0.,0.,0.,f,0.,0.,0.,f);}
inline mat3 dev(const mat3&m){return mat3(m-I3(trace(m)/3.));}
inline mat3 trans(const mat3&m){return mat3(m.xx,m.yx,m.zx,m.xy,m.yy,m.zy,m.xz,m.yz,m.zz);}
inline mat3 inv(const mat3&m){
   const realT det=m.zx*m.xy*m.yz - m.zx*m.xz*m.yy - m.yx*m.xy*m.zz + m.yx*m.xz*m.zy + m.xx*m.yy*m.zz - m.xx*m.yz*m.zy;
   return mat3( ( m.yy*m.zz - m.yz*m.zy) / det,
               -( m.xy*m.zz - m.xz*m.zy) / det,
                ( m.xy*m.yz - m.xz*m.yy) / det,
               -(-m.zx*m.yz + m.yx*m.zz) / det,
                (-m.zx*m.xz + m.xx*m.zz) / det,
               -(-m.yx*m.xz + m.xx*m.yz) / det,
                (-m.zx*m.yy + m.yx*m.zy) / det,
               -(-m.zx*m.xy + m.xx*m.zy) / det,
                (-m.yx*m.xy + m.xx*m.yy) / det);
}

//       stream operations                                 //
inline std::istream& operator>>(std::istream&is,      vec3&v){is>>v.x>>v.y>>v.z;return is;}
inline std::ostream& operator<<(std::ostream&os,const vec3&v){os<<v.x<<'\t'<<v.y<<'\t'<<v.z;return os;}
//inline std::ostream& operator<<(std::ostream&os,const mat3&m){os<<m.xx<<'\t'<<m.xy<<'\t'<<m.xz<<'\t'
                                                                //<<m.yx<<'\t'<<m.yy<<'\t'<<m.yz<<'\t'
                                                                //<<m.zx<<'\t'<<m.zy<<'\t'<<m.zz;return os;}

inline std::ostream& operator<<(std::ostream&os,const mat3&m){
   const realT norm=Frob(m);
   os.precision(6);
   os<<nwl;
   for(indexT i=0;i<3;++i){
      for(indexT j=0;j<3;++j){
         if(abs(m(i,j))/norm>machEps)os<<setw(14)<<m(i,j);
         else os<<setw(14)<<0;
      }
      os<<nwl;
   }
   return os;
}

inline void row6(ostream&os,const mat3&n){
   assert(n.xy==n.yx);
   assert(n.yz==n.zy);
   assert(n.xz==n.zx);
   os.precision(6);
   const int w=14;
   //os<<n.xx<<tab<<n.yy<<tab<<n.zz<<tab<<n.xy<<tab<<n.yz<<tab<<n.xz<<tab;
   os<<setw(w)<<n.xx<<setw(w)<<n.yy<<setw(w)<<n.zz<<setw(w)<<n.xy<<setw(w)<<n.yz<<setw(w)<<n.xz;
}

//       operators involving square root                   //
inline realT mag(const vec3&u){const realT l2=len2(u);return(l2>machEps?std::sqrt(l2):0.);}
inline vec3 unit(const vec3&u){const realT l2=len2(u);assert(l2>machEps);return u/std::sqrt(l2);}
inline realT radius(const vec3&u,const vec3&v){return mag(u-v);}
inline realT vonMises(const mat3&s){
   const realT vm=s.xx*s.xx+s.yy*s.yy+s.zz*s.zz
                  -s.xx*s.yy-s.yy*s.zz-s.zz*s.xx
                  +3.*s.xy*s.yx+3.*s.yz*s.zy+3.*s.zx*s.xz;
   return(vm>machEps?std::sqrt(vm):0.);
}

//       left-to-right operators where order matters       //
inline realT vec3::inner(const vec3&v)const{return x*v.x+y*v.y+z*v.z;}
inline vec3 vec3::inner(const mat3&m)const{return vec3(x*m.xx+y*m.yx+z*m.zx,x*m.xy+y*m.yy+z*m.zy,x*m.xz+y*m.yz+z*m.zz);}
inline mat3 vec3::outer(const vec3&v)const{return mat3(x*v.x,x*v.y,x*v.z,y*v.x,y*v.y,y*v.z,z*v.x,z*v.y,z*v.z);}
inline vec3 mat3::inner(const vec3&v)const{return vec3(xx*v.x+xy*v.y+xz*v.z,yx*v.x+yy*v.y+yz*v.z,zx*v.x+zy*v.y+zz*v.z);}
inline mat3 mat3::inner(const mat3&m)const{return mat3(xx*m.xx+xy*m.yx+xz*m.zx,xx*m.xy+xy*m.yy+xz*m.zy,xx*m.xz+xy*m.yz+xz*m.zz,
                                                       yx*m.xx+yy*m.yx+yz*m.zx,yx*m.xy+yy*m.yy+yz*m.zy,yx*m.xz+yy*m.yz+yz*m.zz,
                                                       zx*m.xx+zy*m.yx+zz*m.zx,xz*m.xy+zy*m.yy+zz*m.zy,zx*m.xz+zy*m.yz+zz*m.zz);}

struct orth6{
   mat3 m;
   vec3 v;
   orth6(){}
   orth6(const mat3&a,const vec3&b):m(a),v(b){}
   orth6&operator=(const orth6&o){m=o.m;v=o.v;return*this;}
   mat3 inner(const mat3&n)const{
      assert(n.xy==n.yx);
      assert(n.yz==n.zy);
      assert(n.xz==n.zx);
      const vec3 s(m.inner(diag(n)));
      const vec3 d(v*vec3(n.xy,n.yz,n.xz));
      return mat3(s.x,d.x,d.z,d.x,s.y,d.y,d.z,d.y,s.z);
   }
};

//ostream&operator<<(ostream&os,const orth6&oo){
   //const double norm=sqrt(dotDot(oo.m,oo.m)+oo.v.inner(oo.v));
   //os.precision(6);
   //os<<nwl;
   //for(int i=0;i<d;++i){
      //for(int j=0;j<d;++j){
         //if(abs(m[i*d+j])/norm>machTol)os<<setw(14)<<m[i*d+j];
         //else os<<setw(14)<<0;
      //}
      //os<<'\n';
   //}
   //return os;
//}


#endif







