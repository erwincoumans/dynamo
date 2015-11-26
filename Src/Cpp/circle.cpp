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
// filename     : circle.cpp
// description	: non-inline methods of class DL_circle
// author	: Bart Barenbrug     June '96
//

#include "circle.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_circle::DL_circle():DL_curve() {
}

void DL_circle::init(DL_geo* geo, DL_point *p, DL_Scalar r, DL_vector *n) {
  DL_curve::init(geo);
  pnt.assign(p);
  // x and y are two vectors each with length r and both perpendicular to
  // each other and n
  if ((n->x!=0) || (n->z!=0)) x.init(n->z,0,-(n->x));
  else {
    if (n->y!=0) x.init(0,n->z,-(n->y));
    else {
      DL_dsystem->get_companion()->Msg("Error: circle::init: zero normal vector\n circle not initialised\n");
      return;
    }
  }
  x.normalize();
  n->crossprod(&x,&y);
  y.normalize();
  x.timesis(r);
  y.timesis(r);
  maxparam=6.283185307;
}

void DL_circle::assign(DL_circle* c) {
  if (c) {
    DL_curve::assign(c);
    pnt.assign(&(c->pnt));
    x.assign(&(c->x));
    y.assign(&(c->y));
  }
}

boolean DL_circle::pos(DL_Scalar s, DL_point *p) {
// returns curve(s) in the point* (in local coordinates of g) and
// whether s is within bounds as return value
  DL_vector v0,v1;
  x.times(cos(s),&v0);
  y.times(sin(s),&v1);
  v0.plusis(&v1);
  pnt.plus(&v0,p);
  return TRUE;
}

boolean DL_circle::deriv(DL_Scalar s, DL_vector *v){
// returns curve'(s) in the vector* (in local coordinates of g)
// and whether s in within bounds as return value
  DL_vector vtmp;
  x.times(-sin(s),v);
  y.times(cos(s),&vtmp);
  v->plusis(&vtmp);
  return TRUE;
}

boolean DL_circle::indomain(DL_Scalar s) {
// returns whether s is within bounds
  return TRUE;
}

DL_Scalar DL_circle::closeto(DL_point *p) {
// returns a curveparameter s with p-curve(s) minimal
  DL_vector vtmp;
  DL_Scalar r,s;
  if (g) {
    DL_point ptemp;
    g->to_local(p,NULL,&ptemp);
    ptemp.minus(&pnt,&vtmp);
  }
  else p->minus(&pnt,&vtmp);
  r=x.norm();
  s=acos(vtmp.inprod(&x)/(r*vtmp.norm()));
  if (vtmp.inprod(&y)<0) s=-s;
  return s;
}
