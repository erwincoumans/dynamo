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
// filename     : torquespring.cpp
// description	: methods of class DL_torquespring
// author	: Bart Barenbrug     April 1998
//

#include "dyna_system.h"
#include "torquespring.h"

// *************** //
// member fuctions //
// *************** //

DL_torquespring::DL_torquespring() {
  c=dc=a=maxtorque=perr=0; g_is_dyna=FALSE;
  d=NULL; g=NULL;
}

DL_torquespring::~DL_torquespring(void) {
}

void DL_torquespring::init(DL_dyna *_d, DL_vector *_dd,
			   DL_vector *_dirv,
                           DL_geo *_g, DL_vector *_dg){
  d=_d; g=_g;
  dd.assign(_dd);
  dg.assign(_dg);
  dirv.assign(_dirv);
  dirv.normalize();
  sens.init(d,&dd,&dirv,g,&dg);
  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;
  storque=0;
  activate();
}

void DL_torquespring::init(DL_dyna *_d, DL_vector *_dd,
                           DL_vector *_dirv,
                           DL_geo *_g, DL_vector *_dg,
			   DL_Scalar _a, DL_Scalar _c) {
  c=_c;
  a=_a;
  init(_d,_dd,_dirv,_g,_dg);
}

void DL_torquespring::init(DL_dyna *_d, DL_vector *_dd,
                           DL_vector *_dirv,
                           DL_geo *_g, DL_vector *_dg,
	                   DL_Scalar _a, DL_Scalar _c, DL_Scalar _dc) {
  dc=_dc;
  init(_d,_dd,_dirv,_g,_dg,_a,_c);
}

void DL_torquespring::calculate_torque(void) {
  DL_Scalar new_error,h=DL_dsystem->get_integrator()->stepsize();

  new_error=sens.sense()-a;

  // update derr
  derr=(new_error-perr)*0.5/h;

  // make perr a little more acurate by estimating its value
  // at t+0.5*h (so you're actually measuring halfway the frame):
  h*=0.5;
  perr=new_error+h*derr;

  // now calculate the actuator value:
  storque=c*perr + dc*derr;
}

void DL_torquespring::calculate_and_apply(void) {
  DL_vector t;
  calculate_torque();
  if (maxtorque!=0) {
    if (maxtorque>0) {
      if (fabs(storque)>maxtorque) {
        deactivate();
        // possibly raise an event here
        return;
      }
    }
    else { // maxforce<0
      DL_Scalar tl=fabs(storque);
      if (tl<maxtorque) {
        storque*=maxtorque/tl;
        // possibly raise an event here
      }
    }
  }
  springtorque(&t);
  d->applytorque(&t);
  if (g_is_dyna) {
    t.neg(&t);
    ((DL_dyna*)g)->applytorque(&t);
  }
}

void DL_torquespring::springtorque(DL_vector *t){
  DL_vector dirvw;
  d->to_world(&dirv,&dirvw);
  dirvw.times(storque,t);
}

void DL_torquespring::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=1;
  if (g_is_dyna) tf=2; else tf=1;
}

void DL_torquespring::get_force_info(int i, DL_actuator_type& at,
			      DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    at=torque;
    _g=d;
    springtorque(_v);
    return;
  }
  if ((i==1) && (g_is_dyna)) {
    at=torque;
    _g=(DL_dyna*)g;
    springtorque(_v);
    _v->neg(_v);
    return;
  }
  at=none;
};
