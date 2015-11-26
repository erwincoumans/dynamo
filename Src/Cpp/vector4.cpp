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

//
// filename: vector4.cpp
// description: quaternions (required by forward dynamics for
//              orientation storage)
//              This class is internal to the DL module
// author: Bart Barenbrug   March '96
//

#include "vector4.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_vector4::DL_vector4(DL_Scalar c0, DL_Scalar c1, DL_Scalar c2, DL_Scalar c3) {
  c[0]=c0; c[1]=c1; c[2]=c2; c[3]=c3;
}

void DL_vector4::init(DL_Scalar n0, DL_Scalar n1, DL_Scalar n2, DL_Scalar n3) {
  c[0]=n0; c[1]=n1; c[2]=n2; c[3]=n3;
}

void DL_vector4::assign(DL_vector4* v) {
  if (v) {
    c[0]=v->c[0]; 
    c[1]=v->c[1]; 
    c[2]=v->c[2]; 
    c[3]=v->c[3];
  }
}

DL_Scalar DL_vector4::norm() {
  return (::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2] + c[3]*c[3]));
}

void DL_vector4::normalize() {
DL_Scalar l=norm();
  if (l!=0.0) {
    l=1.0/l;
    c[0]*=l;
    c[1]*=l;
    c[2]*=l;
    c[3]*=l;
  }
}

DL_Scalar DL_vector4::inprod(DL_vector4* v) {
  if (v) {
	return(c[0]*v->c[0] + c[1]*v->c[1] + c[2]*v->c[2] + c[3]*v->c[3]);
  }
  else return 0;
}

void DL_vector4::plusis(DL_vector4* v) { 
  if (v) {
    c[0] += v->c[0]; 
    c[1] += v->c[1]; 
    c[2] += v->c[2]; 
    c[3] += v->c[3]; 
  }
}

void DL_vector4::timesis(DL_Scalar f) { 
    c[0] *= f; 
    c[1] *= f; 
    c[2] *= f; 
    c[3] *= f; 
}

void DL_vector4::timesis(DL_vector4* v) { 
    c[0] *= v->c[0]; 
    c[1] *= v->c[1]; 
    c[2] *= v->c[2]; 
    c[3] *= v->c[3]; 
}

boolean DL_vector4::equal(DL_vector4* v) {
  if (v) return ((c[0] == v->c[0]) &&
	         (c[1] == v->c[1]) &&
		 (c[2] == v->c[2]) &&
		 (c[3] == v->c[3]));
  else return FALSE;
}

DL_Scalar DL_vector4::get(int r) {
  if ((0<=r) && (r<4)) return c[r];
  else return 0;
}

void DL_vector4::set(int r,DL_Scalar f) {
  if ((0<=r) && (r<4)) c[r]=f;
}


void DL_vector4::plus(DL_vector4* v,DL_vector4 *nv) { 
  nv->c[0]=c[0] + v->c[0];
  nv->c[1]=c[1] + v->c[1];
  nv->c[2]=c[2] + v->c[2];
  nv->c[3]=c[3] + v->c[3];
}

void DL_vector4::minus(DL_vector4 *v,DL_vector4 *nv) { 
  nv->c[0]=c[0] - v->c[0];
  nv->c[1]=c[1] - v->c[1];
  nv->c[2]=c[2] - v->c[2];
  nv->c[3]=c[3] - v->c[3];
}

void DL_vector4::neg(DL_vector4 *nv) {         
  nv->c[0]=-c[0];
  nv->c[1]=-c[1];
  nv->c[2]=-c[2];
  nv->c[3]=-c[3];
}

void DL_vector4::times(DL_Scalar t, DL_vector4 *nv) {
  nv->c[0]=c[0]*t;
  nv->c[1]=c[1]*t;
  nv->c[2]=c[2]*t;
  nv->c[3]=c[3]*t;
}

void DL_vector4::times(DL_vector* w,DL_vector4* nv) {
  //                  1 ( 0  -w^T )
  // implements nv:= ---(         ) self  , where w~ is the outer produkt
  //                  2 ( w  -w~  )         matrix
  //                    (         )

  nv->c[0] = -0.5*(            w->x*c[1] + w->y*c[2] + w->z*c[3]);
  nv->c[1] =  0.5*(w->x*c[0]             + w->z*c[2] - w->y*c[3]);
  nv->c[2] =  0.5*(w->y*c[0] - w->z*c[1]             + w->x*c[3]);
  nv->c[3] =  0.5*(w->z*c[0] + w->y*c[1] - w->x*c[2]            );
}

void DL_vector4::to_matrix(DL_matrix* m){
// PRE: |self|==1
  m->c0.x=2.0*(c[0]*c[0]+c[1]*c[1])-1.0;
    m->c1.x=2.0*(c[1]*c[2]+c[0]*c[3]);
      m->c2.x=2.0*(c[1]*c[3]-c[0]*c[2]);
  m->c0.y=2.0*(c[2]*c[1]-c[0]*c[3]);
    m->c1.y=2.0*(c[0]*c[0]+c[2]*c[2])-1.0;
      m->c2.y=2.0*(c[2]*c[3]+c[0]*c[1]);
  m->c0.z=2.0*(c[3]*c[1]+c[0]*c[2]);
    m->c1.z=2.0*(c[3]*c[2]-c[0]*c[1]);
      m->c2.z=2.0*(c[0]*c[0]+c[3]*c[3])-1.0;
};

void DL_vector4::from_matrix(DL_matrix* m){
 c[0]=::sqrt((m->c0.x+m->c1.y+m->c2.z+1)/4);
 if (c[0]<=0.001) {  // 1st eigenvalue=1, 2nd eigenvalue=-1
   if (m->c0.x > -0.99) {
     c[1]=1+m->c0.x;
     c[2]=m->c1.x;
     c[3]=m->c2.x;
   } else
   if (m->c1.y > -0.99) {
     c[1]=m->c0.y;
     c[2]=1+m->c1.y;
     c[3]=m->c2.y;
   } else
   if (m->c2.z > -0.99) {
     c[1]=m->c0.z;
     c[2]=m->c1.z;
     c[3]=1+m->c2.z;
   }
   normalize();
 }
 else {
   c[1]=(m->c2.y-m->c1.z)/(4*c[0]);
   c[2]=(m->c0.z-m->c2.x)/(4*c[0]);
   c[3]=(m->c1.x-m->c0.y)/(4*c[0]);
 }
};

