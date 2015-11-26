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
// filename: pid.h
// description: a PID controller class
// author: Bart Barenbrug   March '98

#ifndef DL_PIDH
#define DL_PIDH

#include "dyna.h"
#include "controller.h"
#include "sensor.h"
#include "actuator.h"

// ************ //
// class DL_pid //
// ************ //

class DL_pid : public DL_controller {
protected:
  DL_Scalar perr; // the error (for P-term)
  DL_Scalar ierr; // integerated err (for I-term)
  DL_Scalar derr; // derivative of err (for D-term)
  DL_Scalar pcoef,icoef,dcoef; // controller parameters
  DL_Scalar new_error;  // difference between sensor and target
  
  DL_Scalar target; // the controller target

  DL_sensor *sensor;
  DL_actuator* actuator;
public:
  /// for external (to DL) use:
  void	init_coefs(DL_Scalar pc, DL_Scalar ic, DL_Scalar dc);

  void set_target(DL_Scalar t) { target=t; };
  DL_Scalar get_target() { return target; };

  void set_pcoef(DL_Scalar pc){ pcoef=pc; };
  DL_Scalar get_pcoef(){ return pcoef; };
  void set_icoef(DL_Scalar ic){ icoef=ic; };
  DL_Scalar get_icoef(){ return icoef; };
  void set_dcoef(DL_Scalar dc){ dcoef=dc; };
  DL_Scalar get_dcoef(){ return dcoef; };

  // if you want to use an outside actuator and sensor, use
  // the following method to have this controller calculate
  // the actuator value given a provided sensor reading:

  DL_Scalar sens2act(DL_Scalar);
  
             DL_pid(DL_sensor*, DL_actuator*);  // constructor
             DL_pid();
	     ~DL_pid();    // destructor

  virtual void activate(void);     // (re-activate) the controller
  
  /// for (DL) internal use only:
  virtual void calculate_and_apply(void); // calculates and applies the
                                          // controller force
  // for (un)registering with the force_drawer:
  virtual void show_forces(void){ if (actuator) actuator->show_forces(); };
  virtual void hide_forces(void){ if (actuator) actuator->hide_forces(); };

  // for use by the force_drawer:
  virtual void get_fd_info(int& nrf,int& tf){
    if (actuator) actuator->get_fd_info(nrf,tf);
  };
  
  virtual void get_force_info(int i, DL_actuator_type& at,
			      DL_dyna*& d, DL_point *p, DL_vector *v){
    if (actuator) actuator->get_force_info(i,at,d,p,v);
  };
};

#endif
