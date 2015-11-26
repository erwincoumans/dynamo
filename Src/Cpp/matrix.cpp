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
// filename     : matrix.cpp
// description  : non-inline methods of class DL_matrix
// author       : Bart Barenbrug     Januari '97
//
 
#include "matrix.h"
#include "dyna_system.h"
 
// ************************** //
// non-inline member fuctions //
// ************************** //
 
void DL_matrix::makeone() {
  c0.x=c1.y=c2.z=1.0;
  c0.y=c0.z=c1.x=c1.z=c2.x=c2.y=0.0;
}

void DL_matrix::makezero() {
  c0.x=c0.y=c0.z=c1.x=c1.y=c1.z=c2.x=c2.y=c2.z=0.0;
}

void DL_matrix::normalize() {
  c0.normalize();
  c1.normalize();
  c2.normalize();
}

void DL_matrix::assign(DL_matrix* m) {
  c0.assign(&(m->c0));
  c1.assign(&(m->c1));
  c2.assign(&(m->c2));
}

void DL_matrix::assign(DL_vector *v0,DL_vector *v1,DL_vector *v2) {
  c0.assign(v0);
  c1.assign(v1);
  c2.assign(v2);
}

void DL_matrix::assign(DL_Scalar c0x, DL_Scalar c1x, DL_Scalar c2x,
                              DL_Scalar c0y, DL_Scalar c1y, DL_Scalar c2y,
		              DL_Scalar c0z, DL_Scalar c1z, DL_Scalar c2z)
{
  c0.x=c0x; c1.x=c1x; c2.x=c2x;
  c0.y=c0y; c1.y=c1y; c2.y=c2y;
  c0.z=c0z; c1.z=c1z; c2.z=c2z;
}

DL_matrix::DL_matrix(DL_vector* d0, DL_vector* d1, DL_vector* d2) {
 assign(d0,d1,d2);
}

DL_Scalar DL_matrix::get(int r,int c) {
  switch(c) {
    case 0:
      switch(r) {
        case 0: return c0.x;
        case 1: return c0.y;
        case 2: return c0.z;
      };
    case 1:
      switch(r) {
        case 0: return c1.x;
        case 1: return c1.y;
        case 2: return c1.z;
      };
    case 2:
      switch(r) {
        case 0: return c2.x;
        case 1: return c2.y;
        case 2: return c2.z;
      };
  };
  return 0;
}

void DL_matrix::set(int r,int c,DL_Scalar f) {
  switch(c) {
    case 0:
      switch(r) {
        case 0: c0.x=f; break;
        case 1: c0.y=f; break;
        case 2: c0.z=f; break;
      };
      break;
    case 1:
      switch(r) {
        case 0: c1.x=f; break;
        case 1: c1.y=f; break;
        case 2: c1.z=f; break;
      };
      break;
    case 2:
      switch(r) {
        case 0: c2.x=f; break;
        case 1: c2.y=f; break;
        case 2: c2.z=f; break;
      };
      break;
  };
}

void DL_matrix::plus(DL_matrix *m, DL_matrix *nm) {
  c0.plus(&(m->c0), &(nm->c0));
  c1.plus(&(m->c1), &(nm->c1));
  c2.plus(&(m->c2), &(nm->c2));
}

void DL_matrix::minus(DL_matrix *m, DL_matrix *nm) {
  c0.minus(&(m->c0), &(nm->c0));
  c1.minus(&(m->c1), &(nm->c1));
  c2.minus(&(m->c2), &(nm->c2));
}

void DL_matrix::plusis(DL_matrix *m) {
  c0.plusis(&(m->c0));
  c1.plusis(&(m->c1));
  c2.plusis(&(m->c2));
}

void DL_matrix::minusis(DL_matrix *m) {
  c0.minusis(&(m->c0));
  c1.minusis(&(m->c1));
  c2.minusis(&(m->c2));
}

void DL_matrix::timesis(DL_Scalar fac) {
  c0.timesis(fac);
  c1.timesis(fac);
  c2.timesis(fac);
}

void DL_matrix::times(DL_matrix *m, DL_matrix *nm) {
  nm->c0.x=c0.x * m->c0.x + c1.x * m->c0.y + c2.x * m->c0.z;
  nm->c1.x=c0.x * m->c1.x + c1.x * m->c1.y + c2.x * m->c1.z;
  nm->c2.x=c0.x * m->c2.x + c1.x * m->c2.y + c2.x * m->c2.z;

  nm->c0.y=c0.y * m->c0.x + c1.y * m->c0.y + c2.y * m->c0.z;
  nm->c1.y=c0.y * m->c1.x + c1.y * m->c1.y + c2.y * m->c1.z;
  nm->c2.y=c0.y * m->c2.x + c1.y * m->c2.y + c2.y * m->c2.z;

  nm->c0.z=c0.z * m->c0.x + c1.z * m->c0.y + c2.z * m->c0.z;
  nm->c1.z=c0.z * m->c1.x + c1.z * m->c1.y + c2.z * m->c1.z;
  nm->c2.z=c0.z * m->c2.x + c1.z * m->c2.y + c2.z * m->c2.z;
}

void DL_matrix::timestranspose(DL_matrix *m, DL_matrix *nm) {
//nm=self*m^T
  nm->c0.x=c0.x * m->c0.x + c1.x * m->c1.x + c2.x * m->c2.x;
  nm->c1.x=c0.x * m->c0.y + c1.x * m->c1.y + c2.x * m->c2.y;
  nm->c2.x=c0.x * m->c0.z + c1.x * m->c1.z + c2.x * m->c2.z;

  nm->c0.y=c0.y * m->c0.x + c1.y * m->c1.x + c2.y * m->c2.x;
  nm->c1.y=c0.y * m->c0.y + c1.y * m->c1.y + c2.y * m->c2.y;
  nm->c2.y=c0.y * m->c0.z + c1.y * m->c1.z + c2.y * m->c2.z;

  nm->c0.z=c0.z * m->c0.x + c1.z * m->c1.x + c2.z * m->c2.x;
  nm->c1.z=c0.z * m->c0.y + c1.z * m->c1.y + c2.z * m->c2.y;
  nm->c2.z=c0.z * m->c0.z + c1.z * m->c1.z + c2.z * m->c2.z;
}

void DL_matrix::times(DL_Scalar f, DL_matrix *nm) {
  c0.times(f,&(nm->c0));
  c1.times(f,&(nm->c1));
  c2.times(f,&(nm->c2));
}

void DL_matrix::times(DL_vector *v,DL_vector *nv) {
  nv->x=c0.x * v->x + c1.x * v->y + c2.x * v->z;
  nv->y=c0.y * v->x + c1.y * v->y + c2.y * v->z;
  nv->z=c0.z * v->x + c1.z * v->y + c2.z * v->z;
}

void DL_matrix::transposetimes(DL_vector *v,DL_vector *nv) {
  nv->x=c0.x * v->x + c0.y * v->y + c0.z * v->z;
  nv->y=c1.x * v->x + c1.y * v->y + c1.z * v->z;
  nv->z=c2.x * v->x + c2.y * v->y + c2.z * v->z;
}

void DL_matrix::times(DL_point *v,DL_point *nv) {
  nv->x=c0.x * v->x + c1.x * v->y + c2.x * v->z;
  nv->y=c0.y * v->x + c1.y * v->y + c2.y * v->z;
  nv->z=c0.z * v->x + c1.z * v->y + c2.z * v->z;
}

void DL_matrix::tensor(DL_vector *w,DL_vector *v) {
  c0.x=v->x*w->x;    c1.x=v->x*w->y;    c2.x=v->x*w->z;
  c0.y=v->y*w->x;    c1.y=v->y*w->y;    c2.y=v->y*w->z;
  c0.z=v->z*w->x;    c1.z=v->z*w->y;    c2.z=v->z*w->z;
}

void DL_matrix::diag_transpose_vec(DL_vector *d, DL_vector *vec, DL_vector *res){
// res=diag(d)*self^T*vec
  res->init(d->x*(c0.x*vec->x+c0.y*vec->y+c0.z*vec->z),
            d->y*(c1.x*vec->x+c1.y*vec->y+c1.z*vec->z),
            d->z*(c2.x*vec->x+c2.y*vec->y+c2.z*vec->z));
}

void DL_matrix::invert(DL_matrix *nm) {
  DL_Scalar det;
 
  det = c0.x * (c1.y * c2.z - c1.z * c2.y)+
        c1.x * (c2.y * c0.z - c2.z * c0.y)+
        c2.x * (c0.y * c1.z - c0.z * c1.y);
 
  if (det == 0) DL_dsystem->get_companion()->Msg("singular matrix can't be inverted\n");
  else { 
 
    nm->c0.x=(c1.y*c2.z - c1.z*c2.y) / det; 
    nm->c1.x=(c1.z*c2.x - c1.x*c2.z) / det; 
    nm->c2.x=(c1.x*c2.y - c1.y*c2.x) / det;
    nm->c0.y=(c0.z*c2.y - c0.y*c2.z) / det; 
    nm->c1.y=(c0.x*c2.z - c0.z*c2.x) / det; 
    nm->c2.y=(c0.y*c2.x - c0.x*c2.y) / det;
    nm->c0.z=(c0.y*c1.z - c0.z*c1.y) / det; 
    nm->c1.z=(c0.z*c1.x - c0.x*c1.z) / det; 
    nm->c2.z=(c0.x*c1.y - c0.y*c1.x) / det;
  }
}
 
void DL_matrix::transpose(DL_matrix *nm) {
  nm->c0.x=c0.x;  nm->c1.x=c0.y;  nm->c2.x=c0.z;
  nm->c0.y=c1.x;  nm->c1.y=c1.y;  nm->c2.y=c1.z;
  nm->c0.z=c2.x;  nm->c1.z=c2.y;  nm->c2.z=c2.z;
}

void DL_matrix::negcrossdiagcross(DL_point *p, DL_vector *d, DL_point *q){
// self=-crossmatrix(p)*diag(d)*crossmatrix(q)
  c0.x=(p->y*q->y+p->z*q->z)*d->x;
  c0.y=-p->x*q->y*d->y;
  c0.z=-p->x*q->z*d->z;
  c1.x=-p->y*q->x*d->x;
  c1.y=(p->x*q->x+p->z*q->z)*d->y;
  c1.z=-p->y*q->z*d->z;
  c2.x=-p->z*q->x*d->x;
  c2.y=-p->z*q->y*d->y;
  c2.z=(p->x*q->x+p->y*q->y)*d->z;
}

void DL_matrix::diagcrosstranspose(DL_vector *d, DL_point *q, DL_matrix *m){
// self=diag(d)*crossmatrix(q)*m
  c0.x=d->x*(q->y*m->c2.x-q->z*m->c1.x);
  c0.y=d->y*(q->z*m->c0.x-q->x*m->c2.x);
  c0.z=d->z*(q->x*m->c1.x-q->y*m->c0.x);
  c1.x=d->x*(q->y*m->c2.y-q->z*m->c1.y);
  c1.y=d->y*(q->z*m->c0.y-q->x*m->c2.y);
  c1.z=d->z*(q->x*m->c1.y-q->y*m->c0.y);
  c2.x=d->x*(q->y*m->c2.z-q->z*m->c1.z);
  c2.y=d->y*(q->z*m->c0.z-q->x*m->c2.z);
  c2.z=d->z*(q->x*m->c1.z-q->y*m->c0.z);
}

#define ROTATE(a,i,j,k,l) g=a->get(i,j); h=a->get(k,l); \
            a->set(i,j,g-s*(h+g*tau)); a->set(k,l,h+s*(g-h*tau));
 
void DL_matrix::jacobi(DL_matrix* v, DL_vector* dres) {
// jacobi decomposition:
// based on the same function from the source of Walt (waltnumeric2.c).
/* this function decomposes the matrix as 'this = v d vt' 
 * where 'd' is a 3-D vector containing the eigenvalues, 'vt' 
 * is the transpose of 'v', and 'v' is the matrix containing the 
 * eigenvectors (column vectors). The method being used is (of course) 
 * the Jacobi algorithm, for which we refer to standard textbooks. 
 * It uses an iterative approach, iterating over a number of rotations.
 */

 int j,iq,ip,i;
 DL_Scalar tresh, theta, tau, t , sm, s, h, g, c;
 DL_Scalar b[3], d[3], z[3];

 // initialise v to the identity matrix:
 v->makeone();

 // initialise b and d to the diagonal of a and z to 0
 for (ip=0;ip<3;ip++) {
   d[ip]=b[ip]=get(ip,ip);
   z[ip]=0.0;
 }

 // main iteration loop:
 i=0;
 sm=0.0; for (ip=0;ip<3-1;ip++) for (iq=ip+1;iq<3;iq++) sm+=fabs(get(ip,iq));
 while ((i<50) && (sm!=0.0)) {
   if (i<3) tresh=(0.2/9)*sm; else tresh=0.0;
   for (ip=0;ip<2;ip++)
     for (iq=ip+1;iq<3;iq++) {
       g=100.0*fabs(get(ip,iq));
       if (i>3 && fabs(d[ip])+g==fabs(d[ip]) &&
                  fabs(d[iq])+g==fabs(d[iq])  )
         set(ip,iq,0.0);
       else if (fabs(get(ip,iq))>tresh) {
         h=d[iq]-d[ip];
	 if (fabs(h)+g==fabs(h)) t=get(ip,iq)/h;
	 else {
	   theta=0.5*h/get(ip,iq);
	   t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	   if (theta<0.0) t=-t;
	 }
	 c=1.0/sqrt(1.0+t*t);
	 s=t*c;
	 tau=s/(1.0+c);
	 h=t*get(ip,iq); z[ip]-=h; z[iq]+=h;
	 d[ip]-=h; d[iq]+=h; set(ip,iq,0.0);
	 for (j=0;j<=ip-1;j++) {
	   ROTATE(this,j,ip,j,iq)
	 }
	 for (j=ip+1; j<=iq-1; j++) {
	   ROTATE(this,ip,j,j,iq)
	 }
	 for (j=iq+1;j<3;j++) {
	   ROTATE(this,ip,j,iq,j)
	 }
	 for (j=0;j<3;j++) {
	   ROTATE(v,j,ip,j,iq)
	 }
       }
     }
   for (ip=0;ip<3;ip++) {
     b[ip]+=z[ip]; d[ip]=b[ip]; z[ip]=0.0;
   }
   sm=0.0; for (ip=0;ip<3-1;ip++) for (iq=ip+1;iq<3;iq++) sm+=fabs(get(ip,iq));
   i++;
 }
 dres->init(d[0], d[1], d[2]);
}
#undef ROTATE

