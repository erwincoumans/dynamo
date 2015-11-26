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
// filename     : bsplinesegment.cpp
// description	: non-inline methods of class bsplinesegment
// author	: Bart Barenbrug     June '96
//

#include "bsplinesegment.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_bsplinesegment::DL_bsplinesegment():DL_curve() {
  _t=0; _dt=dtinv=1;
}

DL_bsplinesegment::~DL_bsplinesegment(void) {
}

void DL_bsplinesegment::recalc_abcd() {
  a.x=(   - p[0].x + 3.0*p[1].x - 3.0*p[2].x + p[3].x)/6.0;
  a.y=(   - p[0].y + 3.0*p[1].y - 3.0*p[2].y + p[3].y)/6.0;
  a.z=(   - p[0].z + 3.0*p[1].z - 3.0*p[2].z + p[3].z)/6.0;

  b.x=( 3.0*p[0].x - 6.0*p[1].x + 3.0*p[2].x         )/6.0;
  b.y=( 3.0*p[0].y - 6.0*p[1].y + 3.0*p[2].y         )/6.0;
  b.z=( 3.0*p[0].z - 6.0*p[1].z + 3.0*p[2].z         )/6.0;

  c.x=(-3.0*p[0].x +              3.0*p[2].x         )/6.0;
  c.y=(-3.0*p[0].y +              3.0*p[2].y         )/6.0;
  c.z=(-3.0*p[0].z +              3.0*p[2].z         )/6.0;

  d.x=(     p[0].x + 4.0*p[1].x +     p[2].x         )/6.0;
  d.y=(     p[0].y + 4.0*p[1].y +     p[2].y         )/6.0;
  d.z=(     p[0].z + 4.0*p[1].z +     p[2].z         )/6.0;
}

void DL_bsplinesegment::init(DL_geo *geo, DL_point *p0, DL_point *p1, DL_point *p2, DL_point *p3) {
  DL_curve::init(geo);

  p[0].assign(p0); p[1].assign(p1);
  p[2].assign(p2); p[3].assign(p3);

  recalc_abcd();
  
  minparam=_t;
  maxparam=_t+_dt;
}

void DL_bsplinesegment::set_interval(DL_Scalar t, DL_Scalar dt) {
 minparam=_t=t;
 maxparam=t+(_dt=dt);
 dtinv=1.0/dt;
}

void DL_bsplinesegment::assign(DL_bsplinesegment* bss){
  if (bss) {
    DL_curve::assign(bss);
    _t=bss->_t; _dt=bss->_dt; dtinv=bss->dtinv;
    a.assign(&(bss->a));
    b.assign(&(bss->b));
    c.assign(&(bss->c));
    d.assign(&(bss->d));
  }
}

boolean DL_bsplinesegment::pos(DL_Scalar s, DL_point *p) {
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

boolean DL_bsplinesegment::deriv(DL_Scalar s, DL_vector *v){
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

boolean DL_bsplinesegment::deriv2(DL_Scalar s, DL_vector *v){
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

void DL_bsplinesegment::update_control_point(int i, DL_point *pn){
  if ((0<=i) && (i<4)) {
    p[i].assign(pn);
    recalc_abcd();
  }
  else {
    DL_dsystem->get_companion()->Msg("Warning: DL_bsplinesegment::update_control_point(int,DL_point): index out of range\n");
    return;
  }
}

boolean DL_bsplinesegment::indomain(DL_Scalar s) {
// returns whether s is within bounds
  return ((minparam<=s) && (s<=maxparam));
}

DL_Scalar DL_bsplinesegment::closeto(DL_point *p) {
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
