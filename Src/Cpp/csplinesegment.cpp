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
// filename     : csplinesegment.cpp
// description	: non-inline methods of class csplinesegment
// author	: Bart Barenbrug     April '98
//

#include "csplinesegment.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_csplinesegment::DL_csplinesegment():DL_curve() {
  _t=0; _dt=dtinv=1;
}

DL_csplinesegment::~DL_csplinesegment(void) {
}

void DL_csplinesegment::recalc_abcd() {
  float q=0.5;
  a.x=(  -q*p[0].x +(2-q)*p[1].x +  (q-2)*p[2].x +q*p[3].x);
  a.y=(  -q*p[0].y +(2-q)*p[1].y +  (q-2)*p[2].y +q*p[3].y);
  a.z=(  -q*p[0].z +(2-q)*p[1].z +  (q-2)*p[2].z +q*p[3].z);

  b.x=( 2*q*p[0].x +(q-3)*p[1].x +(3-2*q)*p[2].x -q*p[3].x);
  b.y=( 2*q*p[0].y +(q-3)*p[1].y +(3-2*q)*p[2].y -q*p[3].y);
  b.z=( 2*q*p[0].z +(q-3)*p[1].z +(3-2*q)*p[2].z -q*p[3].z);

  c.x=(  -q*p[0].x               +      q*p[2].x          );
  c.y=(  -q*p[0].y               +      q*p[2].y          );
  c.z=(  -q*p[0].z               +      q*p[2].z          );

  d.x=(                   p[1].x                          );
  d.y=(                   p[1].y                          );
  d.z=(                   p[1].z                          );
}

void DL_csplinesegment::init(DL_geo *geo, DL_point *p0, DL_point *p1, DL_point *p2, DL_point *p3) {
  DL_curve::init(geo);

  p[0].assign(p0); p[1].assign(p1);
  p[2].assign(p2); p[3].assign(p3);

  recalc_abcd();
  
  minparam=_t;
  maxparam=_t+_dt;
}

void DL_csplinesegment::set_interval(DL_Scalar t, DL_Scalar dt) {
 minparam=_t=t;
 maxparam=t+(_dt=dt);
 dtinv=1.0/dt;
}

void DL_csplinesegment::assign(DL_csplinesegment* bss){
  if (bss) {
    DL_curve::assign(bss);
    _t=bss->_t; _dt=bss->_dt; dtinv=bss->dtinv;
    a.assign(&(bss->a));
    b.assign(&(bss->b));
    c.assign(&(bss->c));
    d.assign(&(bss->d));
  }
}

boolean DL_csplinesegment::pos(DL_Scalar s, DL_point *p) {
// returns curve(s) in the point* (in local coordinates of g) and
// whether s is within bounds as return value
  DL_Scalar sn=(s-_t)*dtinv; // normalised curve parameter
  boolean inbounds=TRUE;
  if (sn<0) { sn=0; inbounds=FALSE; }
  if (sn>1) { sn=1; inbounds=FALSE; }
  p->x=d.x+sn*(c.x+sn*(b.x+sn*a.x));
  p->y=d.y+sn*(c.y+sn*(b.y+sn*a.y));
  p->z=d.z+sn*(c.z+sn*(b.z+sn*a.z));
  return inbounds;
}

boolean DL_csplinesegment::deriv(DL_Scalar s, DL_vector *v){
// returns curve'(s) in the vector* (in local coordinates of g)
// and whether s in within bounds as return value
  DL_Scalar sn=(s-_t)*dtinv; // normalised curve parameter
  boolean inbounds=TRUE;
  if (sn<0) { sn=0; inbounds=FALSE; }
  if (sn>1) { sn=1; inbounds=FALSE; }
  v->x=c.x+sn*(2*b.x+3*sn*a.x);
  v->y=c.y+sn*(2*b.y+3*sn*a.y);
  v->z=c.z+sn*(2*b.z+3*sn*a.z);
  return inbounds;
}

boolean DL_csplinesegment::deriv2(DL_Scalar s, DL_vector *v){
// returns curve''(s) in the vector* (in local coordinates of g)
// and whether s in within bounds as return value
  DL_Scalar sn=(s-_t)*dtinv; // normalised curve parameter
  boolean inbounds=TRUE;
  if (sn<0) { sn=0; inbounds=FALSE; }
  if (sn>1) { sn=1; inbounds=FALSE; }
  v->x=2*b.x+6*sn*a.x;
  v->y=2*b.y+6*sn*a.y;
  v->z=2*b.z+6*sn*a.z;
  return inbounds;
}

void DL_csplinesegment::update_control_point(int i, DL_point *pn){
  if ((0<=i) && (i<4)) {
    p[i].assign(pn);
    recalc_abcd();
  }
  else {
    DL_dsystem->get_companion()->Msg("Warning: DL_csplinesegment::update_control_point(int,DL_point): index out of range\n");
    return;
  }
}

boolean DL_csplinesegment::indomain(DL_Scalar s) {
// returns whether s is within bounds
  return ((minparam<=s) && (s<=maxparam));
}

DL_Scalar DL_csplinesegment::closeto(DL_point *p) {
// returns a curveparameter s with p-curve(s) minimal
// p given in world coordinates
  // b-splines are rather smooth, so we first take six samples of the
  // curve and then locally use the same algorithm as the ptc::geterror method
  // around the closest of the six samples.
  DL_Scalar s=_t,mins,mind=100000,dist;
  DL_point pl,pc;
  DL_vector pdiff,dc;
  // first convert p to local coordinates
  if (g) g->to_local(p,NULL,&pl);
  else pl.assign(p);
  int i;
  for (i=0; i<=5; i++) {
    pos(s,&pc);
    pc.minus(&pl,&pdiff);
    if ((dist=pdiff.norm())<mind) { mins=s; mind=dist; }
    s+=0.2*_dt;
  }
  s=mins;
  // mins now contains the curve parameter of the six which was closest.
  // now do a local search around s
  #define MAXITER 10
  #define MINDIST (0.00001*_dt)
  i=0;
  do {
    pos(s,&pc); deriv(s,&dc);
    pl.minus(&pc,&pdiff);
    s+=dist=(pdiff.inprod(&dc)/dc.inprod(&dc));
    i++;
  } while ((i<MAXITER) && (dist>MINDIST) && indomain(s));
  
  return s;
}
