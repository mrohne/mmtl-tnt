#ifndef _Vector__h_
#define _Vector__h_

#include <iostream>
#include <iomanip>

class Vector
{
  friend class Matrix;
  friend class Matrix3d;
  friend std::ostream& operator<<(std::ostream&, const Vector&);
  friend std::istream& operator>>(std::istream&, Vector&);
  friend Vector operator*(double, const Vector&);
  friend double max(const Vector&);
  friend double min(const Vector&);
  friend int    max_index(const Vector&);
  friend int    min_index(const Vector&);

public:

  Vector(int = 0);
  Vector(int, double);
  Vector(int, double*);
  Vector(const Vector&);

  ~Vector()    { if(v) delete []v; }

  void free() { if(v) delete []v; }

  int     dim()  const { return d; }
  double* vec()  const { return v; }
  Vector  abs()  const;
  Vector& resize(int);
  void    swap_elem(int,int);

  Vector& operator=(const Vector&);
  double& operator[](int);
  double  operator[](int) const;
  Vector  operator+(const Vector&) const;
  Vector  operator-(const Vector&) const;
  Vector& operator+=(const Vector&);
  Vector& operator-=(const Vector&);
  Vector& operator/=(double);
  Vector& operator*=(double);
  Vector  operator*(double)        const;
  Vector  operator/(double)        const;
  double  operator*(const Vector&) const;
  Vector  operator-() const;

  int     operator==(const Vector&) const;
  int     operator!=(const Vector& w)  const { return !(*this == w); }

private:

  double* v;
  int     d;
};

inline void Print(const Vector& v, std::ostream& out=std::cout) { out << v; }
inline void Read(Vector& v, std::istream& in=std::cin)          { in >> v;  }
//inline double MAX(double x,double y)                  { return (x > y ? x : y); }

double norm(const Vector& v);
double norm2(const Vector& v);
Vector operator^(const Vector& v, const Vector& w);

#endif
