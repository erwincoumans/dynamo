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
// filename: collision.h
// description: collison constraint class
// author: Bart Barenbrug   July '96
//

#ifndef DL_COLLISIONH
#define DL_COLLISIONH

#include "constraint.h"

// ****************** //
// class DL_collision //
// ****************** //

class DL_collision : public DL_constraint {
protected:
  DL_geo *g0,*g1;
  DL_point p0, p1;  // in local coordinates of g0 and g1
  DL_vector v0,v1;  // velocity of contact points at t+h
  DL_point p0w,p1w; // positions (in wc) of contact points at t+h
  DL_Scalar v;          // relative velocity target
  DL_Scalar p;          // relative positional target
  DL_vector n;      // collision normal (in wc)
  boolean g0_is_dyna, g1_is_dyna;
  DL_largematrix dcdX,dXdR;

  void init(DL_geo*, DL_point*,
	    DL_geo*, DL_point*,
	    DL_vector*, int=1);
           // init the constraint with all required parameters

public:
  /// for external use:
             DL_collision(DL_geo*, DL_point*,
	                  DL_geo*, DL_point*,
		          DL_vector*, int=0); // constructor;
	     ~DL_collision(){};        // destructor

  /// for internal (DL) use only:
  // methods for empirical dC/dR determination:
  virtual void begin_test(void);
                     // save the appropriate attributes so they can
                     // be restored after testing, and notify the
		     // dyna's that testing begins.
  virtual void end_test(void);
		     // restore appropriate attributes, and notify
                     // the dyna's that testing has ended

  // methods for analytical dC/dR determination:
  virtual boolean dCdRsub(DL_constraint*, DL_largematrix*);
                     // return the submatrix of dCdR that shows the
		     // relation between the restriction value of this
		     // constraint and the constraint error of
		     // the constraint supplied as parameter
                     // the dimension of sub is cf->dim x dim;
		     // Returns if there is any effect at all (!result=>(m==0))
  virtual boolean dCdFq(DL_dyna*,DL_point*,DL_largematrix*);
                     // calculate the matrix that shows the effect of
		     // application of a (non-central) force to the point of
		     // the dyna on the constraint error of this constraint.
		     // Returns if there is any effect at all (!result=>(m==0))
                     // dcdfq has dimensions dim x 3
  virtual boolean dCdF(DL_dyna*,DL_largematrix*);
                     // calculate the matrix that shows the effect of
		     // application of a central force to the dyna on the
		     // constraint error of this constraint.
		     // Returns if there is any effect at all
		     // (!result=>(m==0))
                     // dcdf has dimensions dim x 3
  virtual boolean dCdM(DL_dyna*,DL_largematrix*);
                     // calculate the matrix that shows the effect of
		     // application of a torque to the dyna on the
		     // constraint error of this constraint.
		     // Returns if there is any effect at all (!result=>(m==0))
                     // dcdm has dimensions dim x 3
  virtual boolean dCdI(DL_dyna*,DL_point*,DL_largematrix*);
                     // calculate the matrix that shows the effect of
		     // application of an impulse to the point of
		     // the dyna on the constraint error of this constraint.
		     // Returns if there is any effect at all (!result=>(m==0))
                     // dcdi has dimensions dim x 3
  
  virtual void apply_restrictions(DL_largevector*);
                     // apply the reaction forces and torques specified
	             // by the parameter
  virtual void get_error(DL_largevector*);
                     // calculate the constraint error vector
  virtual void post_processing(void);
                     // wrap up the calculations for this frame
  // for NOT (un)registering with the force_drawer:
  virtual void show_forces(void){};
  virtual void hide_forces(void){};
};

#endif
