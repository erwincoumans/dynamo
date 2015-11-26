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
// filename     : ellipsoid.cpp
// description	: non-inline methods of class DL_ellipsoid
// author	: Bart Barenbrug     July '96
//

#include "ellipsoid.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_ellipsoid::DL_ellipsoid():DL_surface() {
  minparam0=0;
  maxparam0=6.283185307;
  minparam1=-1.570796327;
  maxparam1=1.570796327;
}

DL_ellipsoid::~DL_ellipsoid() {
}

void DL_ellipsoid::init(DL_geo* geo, DL_point *_c, DL_vector *_x, DL_vector *_y, DL_vector *_z) {
// x,y,z should be orthogonal
  c.assign(_c);
  x.assign(_x);
  y.assign(_y);
  z.assign(_z);
  DL_surface::init(geo);
}

void DL_ellipsoid::assign(DL_ellipsoid* s){
  if (s) {
    DL_surface::assign(s);
    c.assign(&(s->c));
    x.assign(&(s->x));
    y.assign(&(s->y));
    z.assign(&(s->z));
  }
}

boolean DL_ellipsoid::pos(DL_Scalar s, DL_Scalar t, DL_point *p) {
// returns surface(s,t) (in local coordinates of g) in the point* and
// whether s,t is within bounds as return value
  DL_vector v0,v1;
  x.times(sin(s),&v1);
  y.times(cos(s)*cos(t),&v0);
  v1.plus(&v0,&v1);
  z.times(sin(t),&v0);
  v1.plus(&v0,&v1);
  c.plus(&v1,p);
  return TRUE;
}

boolean DL_ellipsoid::deriv0(DL_Scalar s, DL_Scalar t, DL_vector *v){
// returns dsurface(s,t)/ds (in local coordinates of g) in the vector*
// and whether s,t is within bounds as return value
  DL_vector v0;
  x.times(cos(s),v);
  y.times(-sin(s)*cos(t),&v0);
  v->plusis(&v0);
  return TRUE;
}

boolean DL_ellipsoid::deriv1(DL_Scalar s, DL_Scalar t, DL_vector *v){
// returns dsurface(s,t)/ds in the vector* and
// whether s,t is within bounds as return value
  DL_vector v0;
  y.times(-cos(s)*sin(t),v);
  z.times(cos(t),&v0);
  v->plusis(&v0);
  return TRUE;
}

boolean DL_ellipsoid::indomain(DL_Scalar s, DL_Scalar t) {
// returns whether s,t is within bounds
  return TRUE;
}

boolean DL_ellipsoid::closeto(DL_point *p, DL_Scalar& s, DL_Scalar& t) {
// calculates curveparameters s,t with p-surface(s,t) minimal
// returns whether those s,t are within bounds
  DL_vector vtemp, v0;
  if (g) {
    DL_point p0;
    g->to_local(p,NULL,&p0);
    p0.minus(&c,&vtemp);
  }
  else p->minus(&c,&vtemp);
  v0.x=vtemp.inprod(&x)/x.inprod(&x);
  v0.y=vtemp.inprod(&y)/y.inprod(&y);
  v0.z=vtemp.inprod(&z)/z.inprod(&z);
  v0.normalize();
  s=asin(v0.x); if (v0.y<0) s=3.14159-s;
  t=asin(v0.z);
  return TRUE;
}
