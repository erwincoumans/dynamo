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
// filename     : curve.cpp
// description	: non0inline methods of class DL_curve
// author	: Bart Barenbrug     June '96
//

#include "curve.h"
#include "dyna_system.h"
// ************************** //
// non-inline member fuctions //
// ************************** //

DL_curve::DL_curve() {
  g=NULL;
  minparam=-1;
  maxparam=0;
}

void DL_curve::init(DL_geo* geo) {
  g=geo;
}

DL_geo* DL_curve::get_geo(void) {
  return g;
}

void DL_curve::rotatebase(DL_vector *a, DL_vector *b, DL_vector *x, DL_vector *y){
// see constraint::rotatebase
  DL_vector anorm,bnorm,ab,aab,bab,z;
  DL_Scalar za, zab, zaab;
  
  anorm.assign(a); anorm.normalize();
  bnorm.assign(b); bnorm.normalize();
  anorm.crossprod(&bnorm,&ab);
  ab.normalize();
  if (ab.norm()<0.9) return; // no rotation required
  anorm.crossprod(&ab,&aab);
  bnorm.crossprod(&ab,&bab);

  za=x->inprod(&anorm);
  zab=x->inprod(&ab);
  zaab=x->inprod(&aab);
  bnorm.times(za,x);
  ab.times(zab,&z);
  x->plusis(&z);
  bab.times(zaab,&z);
  x->plusis(&z);

  za=y->inprod(&anorm);
  zab=y->inprod(&ab);
  zaab=y->inprod(&aab);
  bnorm.times(za,y);
  ab.times(zab,&z);
  y->plusis(&z);
  bab.times(zaab,&z);
  y->plusis(&z);

  a->assign(b);
}

void DL_curve::assign(DL_curve* c){
  if (c) {
    g=c->g;
    minparam=c->minparam;
    maxparam=c->maxparam;
  }
}

boolean DL_curve::pos(DL_Scalar s, DL_point *p) {
// returns curve(s) in the point* (in local coordinates of g) and
// whether s is within bounds as return value
  p->init(0,0,0);
  DL_dsystem->get_companion()->Msg("Warning: abstract curve::pos(DL_Scalar, point*) called\n");
  return FALSE;
}

boolean DL_curve::deriv(DL_Scalar s, DL_vector *v){
// returns curve'(s) in the vector* (in local coordinates of g)
// and whether s in within bounds as return value
  v->init(0,0,0);
  DL_dsystem->get_companion()->Msg("Warning: abstract curve::deriv(DL_Scalar, vector*) called\n");
  return FALSE;
}

boolean DL_curve::indomain(DL_Scalar s) {
// returns whether s is within bounds
  DL_dsystem->get_companion()->Msg("Warning: abstract curve::indomain(DL_Scalar) called\n");
  return FALSE;
}

DL_Scalar DL_curve::closeto(DL_point *p) {
// returns a curveparameter s with p-curve(s) minimal
  DL_dsystem->get_companion()->Msg("Warning: abstract curve::closeto(point*) called\n");
  return 0;
}
