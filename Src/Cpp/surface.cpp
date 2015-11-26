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
// filename     : surface.cpp
// description	: non-inline methods of class DL_surface
// author	: Bart Barenbrug     June '96
//

#include "surface.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_surface::DL_surface() {
  g=NULL;
  minparam0=minparam1=0;
  maxparam0=maxparam1=0;
}

void DL_surface::init(DL_geo *geo) {
  g=geo;
}

DL_geo* DL_surface::get_geo(void) {
  return g;
}

void DL_surface::assign(DL_surface* s){
  if (s) {
    g=s->g;
    minparam0=s->minparam0;
    minparam1=s->minparam1;
    maxparam0=s->maxparam0;
    maxparam1=s->maxparam1;
  }
}

boolean DL_surface::pos(DL_Scalar s, DL_Scalar t, DL_point *p) {
// returns surface(s,t) (in local coordinates of g) in the point* and
// whether s,t is within bounds as return value
  p->init(0,0,0);
  DL_dsystem->get_companion()->Msg("Warning: abstract surface::pos(DL_Scalar,DL_Scalar,DL_point*) called\n");
  return FALSE;
}

boolean DL_surface::deriv0(DL_Scalar s, DL_Scalar t, DL_vector *v){
// returns dsurface(s,t)/ds (in local coordinates of g) in the vector*
// and whether s,t is within bounds as return value
  v->init(0,0,0);
  DL_dsystem->get_companion()->Msg("Warning: abstract surface::deriv0(DL_Scalar,DL_Scalar,DL_vector*) called\n");
  return FALSE;
}

boolean DL_surface::deriv1(DL_Scalar s, DL_Scalar t, DL_vector *v){
// returns dsurface(s,t)/dt (in local coordinates of g) in the vector*
// and whether s,t is within bounds as return value
  v->init(0,0,0);
  DL_dsystem->get_companion()->Msg("Warning: abstract surface::deriv1(DL_Scalar,DL_Scalar,DL_vector*) called\n");
  return FALSE;
}

boolean DL_surface::indomain(DL_Scalar s, DL_Scalar t) {
// returns whether s,t is within bounds
  DL_dsystem->get_companion()->Msg("Warning: abstract surface::indomain(DL_Scalar,DL_Scalar) called\n");
  return FALSE;
}

boolean DL_surface::closeto(DL_point *p, DL_Scalar& s, DL_Scalar& t) {
// calculates curveparameters s,t with p-surface(s,t) minimal
// returns whether those s,t are within bounds
  DL_dsystem->get_companion()->Msg("Warning: abstract surface::closeto(DL_point*,DL_Scalar&,DL_Scalar&) called\n");
  s=t=0.0;
  return FALSE;
}
