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
// filename: torquespring.h
// description: class describing (damped) torquesprings.
// author: Bart Barenbrug   April 1998

#ifndef DL_TORQUESPRINGH
#define DL_TORQUESPRINGH

#include "dyna.h"
#include "controller.h"
#include "sensor_angle_v.h"

// ********************* //
// class DL_torquespring //
// ********************* //

class DL_torquespring : public DL_controller {
protected:
  DL_dyna *d;  // the torquespring acts between vector dd of d and
  DL_geo *g;   // vector dg of g, around vector dirv
  DL_vector dg,dd,dirv;
  boolean g_is_dyna;
  DL_Scalar perr,derr,storque;
  DL_sensor_angle_v sens;
  void calculate_torque(void); // calculates the spring torque

public:
  /// for external (to DL) use:
  DL_Scalar maxtorque;
               // maximum spring torque before the spring deactivates itself
               // (maxforce: < 0: only apply maximum but do not deactivate
	       //            ==0: no maxtorque checking)
	       //            > 0: deactivate if maximum is exceeded
  DL_Scalar a,c,dc;
               // the rest angle, the spring constant and the damping constant

  void	init(DL_dyna*, DL_vector*, DL_vector*, DL_geo*, DL_vector*);
  void	init(DL_dyna*, DL_vector*, DL_vector*, DL_geo*, DL_vector*,
	     DL_Scalar,DL_Scalar);
  void	init(DL_dyna*, DL_vector*, DL_vector*, DL_geo*, DL_vector*,
	     DL_Scalar,DL_Scalar,DL_Scalar);
        // initialise the spring optionally giving the spring constant,
	// and the rest angle and possibly the damping factor

  void springtorque(DL_vector*); // returns the torquespring force.
     
  virtual void get_fd_info(int&,int&);
  virtual void get_force_info(int, DL_actuator_type&,
		              DL_dyna*&, DL_point*, DL_vector*);

             DL_torquespring();     // constructor
	     ~DL_torquespring();    // destructor

  /// for (DL) internal use only:
  virtual void calculate_and_apply(void); // calculates and applies the
                                          // torquespring force
};

#endif
