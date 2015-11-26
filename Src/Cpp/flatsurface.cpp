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
// filename     : flatsurface.cpp
// description	: non-inline methods of class DL_flatsurface
// author	: Bart Barenbrug     June '96
//

#include "flatsurface.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_flatsurface::DL_flatsurface():DL_surface() {
  minparam0=minparam1=-1000;
  maxparam0=maxparam1=1000;
}

DL_flatsurface::~DL_flatsurface(void) {
}

void DL_flatsurface::init(DL_geo* geo, DL_point *p, DL_vector *n) {
  // dir0 and dir1 are two vectors each with length 1 and both
  // perpendicular to each other and n
  DL_vector d0,d1;
  if ((n->x!=0) || (n->z!=0)) d0.init(n->z,0,-(n->x));
  else {
    if (n->y!=0) d0.init(0,n->z,-(n->y));
    else {
      DL_dsystem->get_companion()->Msg("Error: flatsurface::init: zero normal vector\n flatsurface not initialised\n");
      return;
    }
  }
  n->crossprod(&d0,&d1);
  init(geo,p,&d0,&d1);
}

void DL_flatsurface::init(DL_geo* geo, DL_point *p, DL_vector *d0, DL_vector
*d1) {
  dir0.assign(d0); dir0.normalize();
  dir1.assign(d1); dir1.normalize();
  // check if d0 and d1 are more or less perpendicular:
  if (dir0.inprod(&dir1)>0.1) {
    DL_dsystem->get_companion()->Msg("Warning: flatsurface::init: the two vectors that span the surface should be perpendicular\n");
  }
  pnt.assign(p);
  DL_surface::init(geo);
}

void DL_flatsurface::set_minparam0(DL_Scalar f) {
  minparam0=f;
}

void DL_flatsurface::set_maxparam0(DL_Scalar f) {
  maxparam0=f;
}

void DL_flatsurface::set_minparam1(DL_Scalar f) {
  minparam1=f;
}

void DL_flatsurface::set_maxparam1(DL_Scalar f) {
  maxparam1=f;
}

void DL_flatsurface::assign(DL_flatsurface* s){
  if (s) {
    DL_surface::assign(s);
    pnt.assign(&(s->pnt));
    dir0.assign(&(s->dir0));
    dir1.assign(&(s->dir1));
  }
}

boolean DL_flatsurface::pos(DL_Scalar s, DL_Scalar t, DL_point *p) {
// returns surface(s,t) (in local coordinates of g) in the point* and
// whether s,t is within bounds as return value
  DL_vector v0,v1;
  dir0.times(s,&v0);
  dir1.times(t,&v1);
  v0.plusis(&v1);
  pnt.plus(&v0,p);
  return TRUE;
}

boolean DL_flatsurface::deriv0(DL_Scalar s, DL_Scalar t, DL_vector *v){
// returns dsurface(s,t)/ds (in local coordinates of g) in the vector*
// and whether s,t is within bounds as return value
  v->assign(&dir0);
  return TRUE;
}

boolean DL_flatsurface::deriv1(DL_Scalar s, DL_Scalar t, DL_vector *v){
// returns dsurface(s,t)/ds in the vector* and
// whether s,t is within bounds as return value
  v->assign(&dir1);
  return TRUE;
}

boolean DL_flatsurface::indomain(DL_Scalar s, DL_Scalar t) {
// returns whether s,t is within bounds
  return TRUE;
}

boolean DL_flatsurface::closeto(DL_point *p, DL_Scalar& s, DL_Scalar& t) {
// calculates curveparameters s,t with p-surface(s,t) minimal
// returns whether those s,t are within bounds
  DL_vector vtemp;
  if (g) {
    DL_point p0;
    g->to_local(p,NULL,&p0);
    p0.minus(&pnt,&vtemp);
  }
  else p->minus(&pnt,&vtemp);
  s=dir0.inprod(&vtemp);
  t=dir1.inprod(&vtemp);
  return TRUE;
}
