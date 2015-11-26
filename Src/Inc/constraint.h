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
// filename: constraint.h
// description: class to derive ones own constraints from
// author: Bart Barenbrug   April 1996
//

#ifndef DL_CONSTRAINTH
#define DL_CONSTRAINTH

#include "list.h"
#include "largevector.h"
#include "largematrix.h"
#include "dyna.h"
#include "force_drawable.h"

// ******************* //
// class DL_constraint //
// ******************* //

class DL_constraint : public DL_ListElem, public DL_force_drawable {
protected:
  DL_largevector *F;    // the reaction "force"
  DL_largevector *oldF; // the reaction "force" of the previous frame
  DL_largevector *Fsave;// for saving during testing
  boolean initialised;  // one can only activate an initialised constraint
  boolean testing;   // are we testing?
  boolean veloterms; // add velocityterms to the constraints (to prevent
                     // oscilation in the reaction force). This variable is
		     // set using the detect_osc() method or by the user
		     // thought "soft()" and "hard()"
  boolean veloterms_free; // autodetection for adding veloterms for
                          // this frame
  void detect_osc(); // detect if there are any oscilations and set the
                     // veloterms variable accordingly
  int nr_osc;        // nr frames without veloterms
  int max_osc;       // apply veloterms every one out of max_osc frames
public:
  // methods internal to DL:
  int   dim;         // the dimension of the constraint
  boolean active;    // has the constraint been checked in with constraints
  int	index;       // index used by the constraint manager (the sum of all
                     // dimensions of previous constraint in the list)

  // methods for empirical dC/dR determination:
  virtual void begin_test(void);
                     // save the appropriate attributes so they can
                     // be restored after testing, and notify the
		     // dynas that testing begins.
  virtual void end_test(void);
		     // restore appropriate attributes, and notify
                     // the dynas that testing has ended

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

  virtual void new_frame(void);
                     // prepare for a new frame
  virtual void first_estimate(void);
                     // calculate and apply the first estimate
  virtual void apply_restrictions(DL_largevector*);
                     // apply the reaction forces and torques specified
	             // by the parameter
  virtual boolean check_restrictions();
                     // check the reaction forces and torques
	             // and maybe deactivate
  virtual void test_restriction_changes(DL_largevector*);
                     // announce to the constraint which restriction change
		     // is about to be applied, so the constraint can
		     // decide to deactivate itself
  virtual void apply_restriction_changes(DL_largevector*);
                     // update the restrictions etc. and apply them
  virtual void get_error(DL_largevector*);
                     // calculate the constraint error vector
  virtual void post_processing(void);
                     // wrap up the calculations for this frame

  void reset(void){F->makezero(); oldF->makezero();};
  // reset the current restriction value to 0
  // (to be used after a big global change)
  void reset_undo(void);
  // reset the current restriction value to 0
  // after undoing the current restriction forces
  // (to be used after a big global change)

  void rotatebase(DL_vector*, DL_vector*, DL_vector*, DL_vector*);
                     // rotates the first, third and fourth vectors
		     // in the same way such that the rotated version of
		     // the first vector is the second
		     // Used for rotating local coordinate systems by
		     // for example the ptc and line-hinge constraints

  // methods for external (to DL) use:
  DL_Scalar stiffness;    // the stiffness of the constraint
  void	init(void);       // initialise the constraint
  void  activate(void);   // announce myself at the constraint manager
  void  deactivate(void); // remove myself from the constraint manager
  int   get_dim(void){return dim;} // retrieve the dimension of the constraint
  
  void  soft(void);  // from now on, add velocityterms to soften the constraint
  void  hard(void);  // from now on, do not add velocityterms to harden up the
                     // constraint.
  void  auto_softhard();      // undo the last soft/hard invocation
  void  auto_softhard(int i); // apply velocityterms only every ith frame
  
             DL_constraint();        // constructor
	     ~DL_constraint();       // destructor
};

#endif
