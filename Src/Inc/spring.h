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
// filename: spring.h
// description: class describing (damped) springs.
// author: Bart Barenbrug   August '96

#ifndef DL_SPRINGH
#define DL_SPRINGH

#include "dyna.h"
#include "controller.h"

// *************** //
// class DL_spring //
// *************** //

class DL_spring : public DL_controller {
protected:
  DL_dyna *d;  // the spring attaches point pd of d to point pg of g
  DL_geo *g;
  DL_point pg,pd;
  boolean g_is_dyna;
  DL_vector sforce;

  boolean calculate_force(void); // calculates the spring force

public:
  /// for external (to DL) use:
  DL_Scalar maxforce;
               // maximum spring force before the spring deactivates itself
               // (maxforce: < 0: only appl;y maximum but do not deactivate
	       //            ==0: no maxforce checking)
	       //            > 0: deactivate if maximum is exceeded
  DL_Scalar l,c,dc;
               // the rest length, the spring constant and the damping constant
  boolean el;     // if TRUE: only pulling forces are applied.

  void	init(DL_dyna*, DL_point*, DL_geo*, DL_point*);
  void	init(DL_dyna*, DL_point*, DL_geo*, DL_point*,
	     DL_Scalar,DL_Scalar);
  void	init(DL_dyna*, DL_point*, DL_geo*, DL_point*,
	     DL_Scalar,DL_Scalar,DL_Scalar);
        // initialise the spring optionally giving the spring constant,
	// and the restlength and possibly the damping factor

  void springforce(DL_vector*); // returns the spring force.
     
  virtual void get_fd_info(int&,int&);
  virtual void get_force_info(int, DL_actuator_type&,
		              DL_dyna*&, DL_point*, DL_vector*);

             DL_spring();     // constructor
	     ~DL_spring();    // destructor

  /// for (DL) internal use only:
  virtual void calculate_and_apply(void); // calculates and applies the
                                          // spring force
};

#endif
