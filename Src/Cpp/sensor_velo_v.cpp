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
// filename     : sensor_velo_v.cpp
// description	: methods of class DL_sensor_velo_v
// author	: Bart Barenbrug     March 1998
//

#include "dyna.h"
#include "sensor_velo_v.h"

// *************** //
// member fuctions //
// *************** //

DL_sensor_velo_v::DL_sensor_velo_v() : DL_sensor() {
}

DL_sensor_velo_v::~DL_sensor_velo_v(void) {
}

void DL_sensor_velo_v::init(DL_dyna *_d, DL_point *_pd, DL_vector *_rd,
                            DL_geo  *_g, DL_point* _pg){
  if (_g==_d) {
    DL_dsystem->get_companion()->Msg("Error: sensor_velo_v::init: two different geometries are required!\n sensor not initialised\n");
    return;
  }

  d=_d;
  g=_g;
  pd.assign(_pd);
  pg.assign(_pg);
  rd.assign(_rd);
  rd.normalize();
  if (fabs(rd.norm()-1)>0.001) {
    DL_dsystem->get_companion()->Msg("Error: sensor_velo_v::init: vector should be non-zero\n sensor not initialised\n");
    return;
  }

  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;
}


DL_Scalar DL_sensor_velo_v::sense(){
  g->get_velocity(&pg,&vgw);
  d->get_velocity(&pd,&vdw);
  d->to_world(&rd,&rdw);
  DL_vector pdiff;
  vgw.minus(&vdw,&pdiff);
  return pdiff.inprod(&rdw);
}
