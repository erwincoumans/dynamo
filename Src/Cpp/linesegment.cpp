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
// filename     : linesegment.cpp
// description	: non-inline methods of class linesegment
// author	: Bart Barenbrug     June '96
//

#include "linesegment.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_linesegment::DL_linesegment():DL_line() {
// nothing extra to do...
}

void DL_linesegment::init(DL_geo* geo, DL_point *p0, DL_point *p1) {
  DL_curve::init(geo);
  pnt.assign(p0);
  p1->minus(p0,&dir);
  maxparam=dir.norm();
  dir.normalize();
}

void DL_linesegment::assign(DL_linesegment* s) {
  DL_line::assign(s);
}

boolean DL_linesegment::pos(DL_Scalar s, DL_point *p) {
// returns curve(s) in the point* (in local coordinates of g) and
// whether s is within bounds as return value
  DL_line::pos(s,p);
  return ((0<=s) && (s<=maxparam));
}

boolean DL_linesegment::deriv(DL_Scalar s, DL_vector *v){
// returns curve'(s) in the vector* (in local coordinates of g)
// and whether s in within bounds as return value
  DL_line::deriv(s,v);
  return ((0<=s) && (s<=maxparam));
}

boolean DL_linesegment::indomain(DL_Scalar s) {
// returns whether s is within bounds
  return ((0<=s) && (s<=maxparam));
}

DL_Scalar DL_linesegment::closeto(DL_point *p) {
// returns a curveparameter s with p-curve(s) minimal
  DL_Scalar s=DL_line::closeto(p);
  if (s<0) s=0;
  if (s>maxparam) s=maxparam;
  return s;
}
