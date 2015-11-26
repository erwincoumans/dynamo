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
// filename     : pid.cpp
// description	: methods of class DL_pid
// author	: Bart Barenbrug     Februari 1998
//

#include "dyna_system.h"
#include "pid.h"

// *************** //
// member fuctions //
// *************** //

DL_pid::DL_pid(DL_sensor *s, DL_actuator *a) : DL_controller() {
  if (!s) {
    DL_dsystem->get_companion()->Msg("Error: DL_pid::init: a sensor is required\n PID controller not activated\n");
    sensor=NULL; actuator=NULL;
    return;
  }
  if (!a) {
    DL_dsystem->get_companion()->Msg("Error: DL_pid::init: an actuator is required\n PID controller not activated\n");
    sensor=NULL; actuator=NULL;
    return;
  }

  sensor=s;
  actuator=a;
  perr=sensor->sense();
  ierr=derr=0;
  pcoef=icoef=dcoef=0;
  target=0;

  activate();
}

DL_pid::DL_pid() : DL_controller() {
  perr=ierr=derr=0;
  pcoef=icoef=dcoef=0;
  target=0;
}

DL_pid::~DL_pid(void) {
  deactivate();
}

void DL_pid::activate(void){
  if (sensor && actuator) DL_controller::activate();
  else {
    DL_dsystem->get_companion()->Msg("Error: DL_pid::activate: a sensor and an actuator are required\n PID controller not activated\n");
  }
}

DL_Scalar DL_pid::sens2act(DL_Scalar sens){
  DL_Scalar h=DL_dsystem->get_integrator()->stepsize();

  new_error=sens-target;

  // update derr
  derr=(new_error-perr)*0.5/h;

  // update ierr:
  ierr+=new_error;

  // make perr a little more acurate by estimating its value
  // at t+0.5*h (so you're actually measuring halfway the frame):
  h*=0.5;
  perr=new_error+h*derr;

  // now calculate and apply the actuator values:
  return perr*pcoef + ierr*icoef + derr*dcoef;
}

void DL_pid::calculate_and_apply(void) {
  actuator->apply(sens2act(sensor->sense()));
}

void DL_pid::init_coefs(DL_Scalar pc, DL_Scalar ic, DL_Scalar dc){
  pcoef=pc;
  icoef=ic;
  dcoef=dc;
}
