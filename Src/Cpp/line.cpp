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
// filename     : line.cpp
// description	: non-inline methods of class DL_line
// author	: Bart Barenbrug     June '96
//

#include "line.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_line::DL_line():DL_curve() {
}

DL_line::~DL_line() {
}

void DL_line::init(DL_geo* geo, DL_point *p, DL_vector *d) {
  DL_curve::init(geo);
  pnt.assign(p);
  dir.assign(d);
  dir.normalize();
  minparam=-(maxparam=1000);
}

void DL_line::set_minparam(DL_Scalar f) {
  minparam=f;
}

void DL_line::set_maxparam(DL_Scalar f) {
  maxparam=f;
}
  
void DL_line::assign(DL_line* s) {
  if (s) {
    DL_curve::assign(s);
    pnt.assign(&(s->pnt));
    dir.assign(&(s->dir));
  }
}

boolean DL_line::pos(DL_Scalar s, DL_point *p) {
// returns curve(s) in the point* (in local coordinates of g) and
// whether s is within bounds as return value
  DL_vector v;
  dir.times(s,&v);
  pnt.plus(&v,p);
  return TRUE;
}

boolean DL_line::deriv(DL_Scalar s, DL_vector *v){
// returns curve'(s) in the vector* (in local coordinates of g)
// and whether s in within bounds as return value
  v->assign(&dir);
  return TRUE;
}

boolean DL_line::indomain(DL_Scalar s) {
// returns whether s is within bounds
  return TRUE;
}

DL_Scalar DL_line::closeto(DL_point *p) {
// returns a curveparameter s with p-curve(s) minimal
  DL_vector vtemp;
  if (g) {
    DL_point p0;
    g->to_local(p,NULL,&p0);
    p0.minus(&pnt,&vtemp);
  }
  else p->minus(&pnt,&vtemp);
  return dir.inprod(&vtemp);
}
