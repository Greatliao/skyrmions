// By J. Burgy Jan 30 2001
#ifndef COMPLEX_H_
#define COMPLEX_H_

#include <math.h>
#include <iostream.h>

class complex
{
  //friend long double abs(const complex &);
  //friend long double norm(const complex &);
  //friend long double arg(const complex &);
  //friend complex conj(const complex &);
  // friend complex exp(const complex &);
  //friend long double imag(const complex &);
  //friend long double real(const complex &);
  // Modification by G. Alvarez
  //friend complex operator*(double, const complex &);
  // end of modification
  //friend ostream & operator<<(ostream &, complex &);
  //friend istream & operator>>(istream &, complex &);


 
 public:
  
  long double re;
  long double im;
  inline complex(long double, long double = 0.0);
  inline complex();

  complex operator+(const complex &);
  complex operator-(const complex &);
  complex operator*(const complex &);
  complex operator*(const long double &);
  complex operator/(const long double &);
  complex & operator+=(const complex &);
  complex & operator-=(const complex &);
  complex & operator*=(const complex &);
  complex & operator/=(const long double &);
  bool operator==(const complex &);
  bool operator!=(const complex &);
};

static const complex complex_zero(0.0,0.0);

inline complex::complex(long double x, long double y) : re(x), im(y) {}
inline complex::complex () : re(0.0), im(0.0) {}
inline long double norm(const complex &z){ return z.re*z.re+z.im*z.im; }
inline long double abs(const complex &z){ return (norm(z)); }
inline complex conj(const complex &z){ return complex(z.re,-z.im); }
inline long double imag(const complex &z){ return z.im; }
inline long double real(const complex &z){ return z.re; }
inline complex complex::operator+(const complex &z){ 
  return complex(re+z.re,im+z.im); }
inline complex complex::operator-(const complex &z){ 
  return complex(re-z.re,im-z.im); }
inline complex complex::operator*(const complex &z){
  return complex(re*z.re-im*z.im,re*z.im+im*z.re); }
inline complex complex::operator*(const long double &d){
  return complex(re*d,im*d); }
inline complex complex::operator/(const long double &d){
  return complex(re/d,im/d); }
inline complex & complex::operator+=(const complex &z){
  re += z.re; im += z.im; return *this; }
inline complex & complex::operator-=(const complex &z){
  re -= z.re; im -= z.im; return *this; }
inline complex & complex::operator/=(const long double &d){
  re /= d; im /= d; return *this; }
inline bool complex::operator==(const complex &z){ 
  return (re==z.re && im == z.im); }
inline bool complex::operator!=(const complex &z){ 
  return (re!=z.re || im != z.im); }


inline istream & operator>>(istream &is, complex &z){
  long double x = 0.0, y = 0.0;
  char c;

  is >> c;
  if(c == '('){
    is >> x >> c;
    if(c == ',') { is >> y >> c; }
    if(c != ')') { is.clear(ios::failbit); }
  } else {
    is.putback(c);
    is >> x;
  }

  z = complex(x,y);
  return(is);
}

// Modification by G. Alvarez
complex operator*(double d, const complex &c);
complex exp(const complex &c);

#endif
