// Philip Wallstedt 2004-2009
#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <iomanip>
#include <cassert>
#include <limits>
#include <cmath>
const double pi = 3.1415926535897932384626433832795028841971694;
const double machTol = 16.0 * std::numeric_limits<double>::epsilon();
using namespace std;

template<class charT, class traits>
basic_ostream<charT, traits>& tab(basic_ostream<charT, traits>& os) {os << '\t' << setw(10); return os;}
template<class charT, class traits>
basic_ostream<charT, traits>& nwl(basic_ostream<charT, traits>& os) {os << '\n' << setw(10); return os;}

// The gcc optimization -O3 gives bad results for ceil unless
// the machTol is used in this way.
inline double Pceil(const double& d) {return ceil(d - machTol);}

template<typename T>inline T abs(T a) {return (a < 0. ? -a : a);}
template<typename T>inline T sgn(T a) {return (a < 0. ? -1. : 1.);}
template<bool>struct Assert;
template<>struct Assert<true> {}; // Assert<(1==2)>();

template<typename T>T max(const vector<T>& v) {
   double m = v[0]; for (unsigned i = 0; i < v.size(); i += 1)m = (v[i] > m ? v[i] : m); return m;
}

template<typename T>T min(const vector<T>& v) {
   double m = v[0]; for (unsigned i = 0; i < v.size(); i += 1)m = (v[i] < m ? v[i] : m); return m;
}

struct Matrix2;

struct Vector2 {
   double x, y;
   Vector2() {x = 0.; y = 0.;}
   Vector2(double a, double b) {x = a; y = b;}
   Vector2(const Vector2& v) {x = v.x; y = v.y;}
   Vector2& operator=(const double& a) {x = a; y = a; return *this;}
   Vector2& operator=(const Vector2& v) {x = v.x; y = v.y; return *this;}
   Vector2& operator+=(const Vector2& v) {x += v.x; y += v.y; return *this;}
   Vector2& operator-=(const Vector2& v) {x -= v.x; y -= v.y; return *this;}
   Vector2& operator*=(const Vector2& v) {x *= v.x; y *= v.y; return *this;}
   Vector2& operator/=(const Vector2& v) {x /= v.x; y /= v.y; return *this;}
   Vector2& operator+=(const double& v) {x += v; y += v; return *this;}
   Vector2& operator-=(const double& v) {x -= v; y -= v; return *this;}
   Vector2& operator*=(const double& v) {x *= v; y *= v; return *this;}
   Vector2& operator/=(const double& v) {x /= v; y /= v; return *this;}
   const Vector2 operator-(               )const {return Vector2(-x, -y);}
   const Vector2 operator+(const Vector2& v)const {return Vector2(x + v.x, y + v.y);}
   const Vector2 operator-(const Vector2& v)const {return Vector2(x - v.x, y - v.y);}
   const Vector2 operator*(const Vector2& v)const {return Vector2(x * v.x, y * v.y);}
   const Vector2 operator/(const Vector2& v)const {return Vector2(x / v.x, y / v.y);}
   const Vector2 operator+(const double& a)const {return Vector2(x + a, y + a);}
   const Vector2 operator-(const double& a)const {return Vector2(x - a, y - a);}
   const Vector2 operator*(const double& a)const {return Vector2(x * a, y * a);}
   const Vector2 operator/(const double& a)const {return Vector2(x / a, y / a);}
   bool operator==(const Vector2& v) {return x == v.x && y == v.y;}
   bool operator==(const double & v) {return x == v  && y == v  ;}
   bool operator!=(const Vector2& v) {return x != v.x || y != v.y;}
   bool operator!=(const double & v) {return x != v  || y != v  ;}
   inline double  inner(const Vector2& v)const;
   inline Vector2 inner(const Matrix2& m)const;
   inline Matrix2 outer(const Vector2& v)const;
};

inline const Vector2 operator+(const double& a, const Vector2& v) {return Vector2(a + v.x, a + v.y);}
inline const Vector2 operator-(const double& a, const Vector2& v) {return Vector2(a - v.x, a - v.y);}
inline const Vector2 operator*(const double& a, const Vector2& v) {return Vector2(a * v.x, a * v.y);}
inline const Vector2 operator/(const double& a, const Vector2& v) {return Vector2(a / v.x, a / v.y);}

struct Matrix2 {
   double xx, xy, yx, yy;
   Matrix2() {xx = 0.; xy = 0.; yx = 0.; yy = 0.;}
   Matrix2(double a) {xx = a; xy = a; yx = a; yy = a;}
   Matrix2(double a, double b, double c, double d) {xx = a; xy = b; yx = c; yy = d;}
   Matrix2(const Matrix2& m) {xx = m.xx; xy = m.xy; yx = m.yx; yy = m.yy;}
   Matrix2& operator=(const double& a) {xx = a; xy = a; yx = a; yy = a; return *this;}
   Matrix2& operator=(const Matrix2& m) {xx = m.xx; xy = m.xy; yx = m.yx; yy = m.yy; return *this;}
   Matrix2& operator+=(const Matrix2& m) {xx += m.xx; xy += m.xy; yx += m.yx; yy += m.yy; return *this;}
   Matrix2& operator-=(const Matrix2& m) {xx -= m.xx; xy -= m.xy; yx -= m.yx; yy -= m.yy; return *this;}
   Matrix2& operator*=(const Matrix2& m) {xx *= m.xx; xy *= m.xy; yx *= m.yx; yy *= m.yy; return *this;}
   Matrix2& operator/=(const Matrix2& m) {xx /= m.xx; xy /= m.xy; yx /= m.yx; yy /= m.yy; return *this;}
   Matrix2& operator+=(const double& m) {xx += m; xy += m; yx += m; yy += m; return *this;}
   Matrix2& operator-=(const double& m) {xx -= m; xy -= m; yx -= m; yy -= m; return *this;}
   Matrix2& operator*=(const double& m) {xx *= m; xy *= m; yx *= m; yy *= m; return *this;}
   Matrix2& operator/=(const double& m) {xx /= m; xy /= m; yx /= m; yy /= m; return *this;}
   const Matrix2 operator+(const Matrix2& m)const {return Matrix2(xx + m.xx, xy + m.xy, yx + m.yx, yy + m.yy);}
   const Matrix2 operator-(const Matrix2& m)const {return Matrix2(xx - m.xx, xy - m.xy, yx - m.yx, yy - m.yy);}
   const Matrix2 operator*(const Matrix2& m)const {return Matrix2(xx * m.xx, xy * m.xy, yx * m.yx, yy * m.yy);}
   const Matrix2 operator/(const Matrix2& m)const {return Matrix2(xx / m.xx, xy / m.xy, yx / m.yx, yy / m.yy);}
   const Matrix2 operator+(const double& a)const {return Matrix2(xx + a, xy + a, yx + a, yy + a);}
   const Matrix2 operator-(const double& a)const {return Matrix2(xx - a, xy - a, yx - a, yy - a);}
   const Matrix2 operator*(const double& a)const {return Matrix2(xx * a, xy * a, yx * a, yy * a);}
   const Matrix2 operator/(const double& a)const {return Matrix2(xx / a, xy / a, yx / a, yy / a);}
   bool operator==(const Matrix2& m) {return xx == m.xx && xy == m.xy && yx == m.yx && yy == m.yy;}
   bool operator!=(const Matrix2& m) {return xx != m.xx || xy != m.xy || yx != m.yx || yy != m.yy;}
   inline Vector2 inner(const Vector2& v)const;
   inline Matrix2 inner(const Matrix2& n)const;
};
inline const Matrix2 operator+(const double& a, const Matrix2& m) {return Matrix2(a + m.xx, a + m.xy, a + m.yx, a + m.yy);}
inline const Matrix2 operator-(const double& a, const Matrix2& m) {return Matrix2(a - m.xx, a - m.xy, a - m.yx, a - m.yy);}
inline const Matrix2 operator*(const double& a, const Matrix2& m) {return Matrix2(a * m.xx, a * m.xy, a * m.yx, a * m.yy);}
inline const Matrix2 operator/(const double& a, const Matrix2& m) {return Matrix2(a / m.xx, a / m.xy, a / m.yx, a / m.yy);}

//       unary operators                                   //
inline double len2(const Vector2& u) {return u.x * u.x + u.y * u.y;}
inline double trace(const Matrix2& m) {return m.xx + m.yy;}
inline double det(const Matrix2& m) {return m.xx * m.yy - m.xy * m.yx;}
inline double absMax(const Matrix2& s) {return max(max(abs(s.xx), abs(s.xy)), max(abs(s.yx), abs(s.yy)));}
inline Matrix2 I2(const double f = 1.0) {return Matrix2(f, 0.0, 0.0, f);}
inline Matrix2 trans(const Matrix2& m) {return Matrix2(m.xx, m.yx, m.xy, m.yy);}
inline Matrix2 inv(const Matrix2& m) {const double d = det(m); return Matrix2(m.yy / d, -m.xy / d, -m.yx / d, m.xx / d);}

//       printing operations                               //
inline std::ostream& operator<<(std::ostream& os, const Vector2& v) {os << v.x << tab << v.y; return os;}
inline std::istream& operator>>(std::istream& is,      Vector2& v) {is >> v.x     >> v.y; return is;}
inline std::ostream& operator<<(std::ostream& os, const Matrix2& m) {os << m.xx << tab << m.xy << tab << m.yx << tab << m.yy; return os;}
inline std::istream& operator>>(std::istream& is,      Matrix2& v) {is >> v.xx     >> v.xy     >> v.yx     >> v.yy; return is;}

//       operators involving square root                   //
inline double mag(const Vector2& u) {return std::sqrt(u.x * u.x + u.y * u.y);}
inline Vector2 unit(const Vector2& u) {assert(abs(u.x) > machTol || abs(u.y) > machTol); return u / mag(u);}
inline double radius(const Vector2& u, const Vector2& v) {
   const double s = (u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y);
   return (s > machTol ? std::sqrt(s) : 0.) ;
}
inline double vonMises(const Matrix2& s) {return sqrt(s.xx * s.xx - s.xx * s.yy + s.yy * s.yy + 3.*s.xy * s.yx);}

//       left-to-right operators where order matters       //
inline double  Vector2::inner(const Vector2& v)const {return x * v.x + y * v.y;}
inline Vector2 Vector2::inner(const Matrix2& m)const {return Vector2(m.xx * x + m.yx * y, m.xy * x + m.yy * y);}
inline Matrix2 Vector2::outer(const Vector2& v)const {return Matrix2(x * v.x, x * v.y, y * v.x, y * v.y);}
inline Vector2 Matrix2::inner(const Vector2& v)const {return Vector2(xx * v.x + xy * v.y, yx * v.x + yy * v.y);}
inline Matrix2 Matrix2::inner(const Matrix2& n)const {
   return Matrix2(xx * n.xx + xy * n.yx, xx * n.xy + xy * n.yy,
                  yx * n.xx + yy * n.yx, yx * n.xy + yy * n.yy);
}

#endif