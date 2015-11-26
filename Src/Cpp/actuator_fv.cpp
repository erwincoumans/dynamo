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
// filename     : actuator_fv.cpp
// description	: methods of class DL_actuator_fv
// author	: Bart Barenbrug     Februari 1998
//

#include "dyna.h"
#include "actuator_fv.h"

// *************** //
// member fuctions //
// *************** //

DL_actuator_fv::DL_actuator_fv() : DL_actuator() {
}

DL_actuator_fv::~DL_actuator_fv(void) {
}

void DL_actuator_fv::init(DL_dyna *_d, DL_point *_pd, DL_vector *_rd,
                          DL_geo  *_g, DL_point* _pg){
  if (_g==_d) {
    DL_dsystem->get_companion()->Msg("Error: actuator_fv::init: two different geometries are required!\n actuator not initialised\n");
    return;
  }

  d=_d;
  g=_g;
  pd.assign(_pd);
  pg.assign(_pg);
  rd.assign(_rd);
  rd.normalize();
  if (fabs(rd.norm()-1)>0.001) {
    DL_dsystem->get_companion()->Msg("Error: actuator_fv::init: vector should be non-zero\n actuator not initialised\n");
    return;
  }

  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;
  
}


void DL_actuator_fv::apply(DL_Scalar a){
  d->to_world(&rd,&f);
  f.timesis(limit(a));
  d->applyforce(&pd,d,&f);
  if (g_is_dyna) {
    DL_vector mf;
    f.neg(&mf);
    ((DL_dyna*)g)->applyforce(&pg,g,&mf);
  }
}

void DL_actuator_fv::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=1;
  if (g_is_dyna) tf=2; else tf=1;
}

void DL_actuator_fv::get_force_info(int i, DL_actuator_type& at,
			            DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    at=force;
    _g=d;
    _p->assign(&pd);
    _v->assign(&f);
    return;
  }
  if ((i==1) && (g_is_dyna)) {
    at=force;
    _g=(DL_dyna*)g;
    _p->assign(&pg);
    f.neg(_v);
    return;
  }
  at=none;
};
