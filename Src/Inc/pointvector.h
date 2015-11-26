/*
  DYNAMO - Dynamic Motion library
  Copyright (C) 1996-1999 Bart Barenbrug

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Library General Public
  License as published by the Free Software Foundation; either
  version 2 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Library General Public License for more details.

  You should have received a copy of the GNU Library General Public
  License along with this library; if not, write to the Free
  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  Please send remarks, questions and bug reports to bartb@win.tue.nl,
  or write to:
                  Bart Barenbrug
		  Department of Mathematics and Computing Science
		  Eindhoven University of Technology
		  P.O. Box 513, 5600 MB Eindhoven, The Netherlands
*/

// filename		: pointvector.h
// description: generic point and vector classes
// author: Bart Barenbrug   March '96
//

#ifndef DL_POINTVECTORH 
#define DL_POINTVECTORH

#include <iostream.h>
#include <math.h>
#include "boolean.h"
#include "scalar.h"
#include "list.h"

class DL_point;
class DL_vector;

// *************** //
// class DL_vector //
// *************** //

class DL_vector {
  public:   
    DL_Scalar  x,y,z;

    void       init(DL_Scalar,DL_Scalar,DL_Scalar);
    void       assign(DL_vector*);
    DL_Scalar  norm();
    void       normalize();
    DL_Scalar  inprod(DL_vector*);
    void       plusis(DL_vector*);
    void       minusis(DL_vector*);
    void       timesis(DL_Scalar);
    boolean    equal(DL_vector*);

    DL_Scalar  get(int);
    void       set(int,DL_Scalar);

    void       crossprod(DL_vector*,DL_vector*);   
    void       neg(DL_vector*);
    void       times(DL_Scalar,DL_vector*);       
    void       plus(DL_vector*,DL_vector*);      
    void       minus(DL_vector*,DL_vector*);     
    void       topoint(DL_point*);

               DL_vector();                // constructor
               DL_vector(DL_vector*);      // copy constructor
               DL_vector(DL_Scalar,DL_Scalar,DL_Scalar);  // constructor
               ~DL_vector(){};                // destructor
};

// ************** //
// class DL_point //
// ************** //

class DL_point : public DL_ListElem {
  public:   
    DL_Scalar  x,y,z;

    void       init(DL_Scalar,DL_Scalar,DL_Scalar);
    void       assign(DL_point*);
    boolean    equal(DL_point*);

    void       plus(DL_vector*,DL_point*);
    void       minus(DL_point*,DL_vector*);
    void       times(DL_Scalar,DL_point*);
    void       plusis(DL_vector*);
    void       minusis(DL_vector*);
    void       timesis(DL_Scalar);
    void       tovector(DL_vector*);

               DL_point();                  // constructor
	       DL_point(DL_point*);         // copy constructor
               DL_point(DL_Scalar,DL_Scalar,DL_Scalar); // constructor
	       ~DL_point(){};               // destructor
};

inline DL_point::DL_point():DL_ListElem() {
//  x=y=z=0;
};

inline DL_point::DL_point(DL_Scalar nx, DL_Scalar ny, DL_Scalar nz) {
  x=nx; y=ny; z=nz;
};

inline void DL_point::init(DL_Scalar nx, DL_Scalar ny, DL_Scalar nz) {
  x=nx; y=ny; z=nz;
};

inline void DL_point::assign(DL_point* p) {
  if (p) {
    x=p->x;
    y=p->y;
    z=p->z;
  }
  else { 
    x=y=z=0; 
  }
};

inline DL_point::DL_point(DL_point* p) {
  assign(p);
};

inline void DL_point::plus(DL_vector* v,DL_point* nv) { 
  nv->x=x + v->x;
  nv->y=y + v->y;
  nv->z=z + v->z;
};

inline void DL_point::minus(DL_point* v,DL_vector* nv) { 
  nv->x=x - v->x;
  nv->y=y - v->y;
  nv->z=z - v->z;
};

inline void DL_point::times(DL_Scalar f,DL_point* nv) { 
  nv->x=x * f;
  nv->y=y * f;
  nv->z=z * f;
};

inline void DL_point::plusis(DL_vector* v) { 
  x+=v->x;
  y+=v->y;
  z+=v->z;
};

inline void DL_point::minusis(DL_vector* v) { 
  x-=v->x;
  y-=v->y;
  z-=v->z;
};

inline void DL_point::timesis(DL_Scalar f) { 
  x*=f;
  y*=f;
  z*=f;
};

inline boolean DL_point::equal(DL_point* v) {
  return ( (x == v->x) && (y == v->y) && (z == v->z));
};

inline void DL_point::tovector(DL_vector *nv) {
  nv->x=x;
  nv->y=y;
  nv->z=z;
};


inline DL_vector::DL_vector() {
//  x=y=z=0;
};

inline void DL_vector::assign(DL_vector* v) {
  x=v->x; y=v->y; z=v->z;
};

inline DL_vector::DL_vector(DL_vector* v) {
  assign(v);
};

inline DL_vector::DL_vector(DL_Scalar nx, DL_Scalar ny, DL_Scalar nz) {
  x=nx; y=ny; z=nz;
};

inline void DL_vector::init(DL_Scalar nx, DL_Scalar ny, DL_Scalar nz) {
  x=nx; y=ny; z=nz;
};


inline DL_Scalar DL_vector::norm() {
  return (::sqrt(x*x + y*y + z*z));
};

inline void DL_vector::normalize() {
  DL_Scalar l=norm();
  if (l!=0.0) {
    x/=l;
    y/=l;
    z/=l;
  }
};

inline void DL_vector::plusis(DL_vector* v) { 
    x = x+v->x; 
    y = y+v->y; 
    z = z+v->z;
};

inline void DL_vector::minusis(DL_vector* v) { 
    x = x-v->x; 
    y = y-v->y; 
    z = z-v->z;
};

inline void DL_vector::timesis(DL_Scalar f) { 
    x = x*f; 
    y = y*f; 
    z = z*f;
};

inline boolean DL_vector::equal(DL_vector* v) {
  return ( (x == v->x) && (y == v->y) && (z == v->z));
};

inline DL_Scalar DL_vector::get(int r) {
  switch(r) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
  }
  return 0;
};

inline void DL_vector::set(int r,DL_Scalar f) {
  switch(r) {
    case 0: x=f;
            break;
    case 1: y=f;
	    break;
    case 2: z=f;
	    break;
  }
};


inline void DL_vector::plus(DL_vector* v,DL_vector *nv) { 
  nv->x=x + v->x;
  nv->y=y + v->y; 
  nv->z=z + v->z;
};

inline void DL_vector::minus(DL_vector *v,DL_vector *nv) { 
  nv->x=x - v->x; 
  nv->y=y - v->y; 
  nv->z=z - v->z;
};

inline void DL_vector::neg(DL_vector *nv) {         
  nv->x=-x;
  nv->y=-y;
  nv->z=-z;
}

inline void DL_vector::times(DL_Scalar t, DL_vector *nv) {
  nv->x=x*t;
  nv->y=y*t;
  nv->z=z*t;
};

inline void DL_vector::crossprod(DL_vector *v,DL_vector *nv) {
  nv->x=y*v->z - z*v->y;
  nv->y=z*v->x - x*v->z;
  nv->z=x*v->y - y*v->x;
};

inline DL_Scalar DL_vector::inprod(DL_vector *v) {
  return x*v->x + y*v->y + z*v->z;
};


inline void DL_vector::topoint(DL_point *np) {
  np->x=x;
  np->y=y;
  np->z=z;
};

#endif
