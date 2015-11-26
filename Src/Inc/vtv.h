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
// filename: vtv.h
// description: velocity-of-point to velocity-of-point constraint class
// author: Bart Barenbrug   July '96
//

#ifndef DL_VTVH
#define DL_VTVH

#include "constraint.h"
#include "constraint_manager.h"

// ************ //
// class DL_vtv //
// ************ //

class DL_vtv : public DL_constraint {
protected:
  DL_dyna  *d;   // the constraint connects point pd in d to
  DL_point pd;  // point pg in g
  DL_geo   *g;
  DL_point pg;
  DL_point pgw; // pg in wc (only calculated once if !g_is_dyna)
  DL_vector dpgw;  // velocity of pg (in wc; used only if !g_is_dyna) 
  boolean g_is_dyna;
  
public:
  // for external access:
  DL_Scalar maxforce; // maximum force (<=0: no maxforce checking)
  void init(DL_dyna*, DL_point*, DL_geo*, DL_point*);
           // initialise the constraint
  void reactionforce(DL_vector*);
           // returns a copy of the reaction force (in world coordinates)

  DL_dyna* get_dyna(void){ return d; };
  void get_dyna_point(DL_point *dp){ dp->assign(&pd); };
  void set_dyna_point(DL_point *dp){ pd.assign(dp); };
  DL_geo* get_geo(void){ return g; };
  void get_geo_point(DL_point *gp){ gp->assign(&pg); };
  void set_geo_point(DL_point *gp){ pg.assign(gp); };

  virtual void get_fd_info(int&,int&);
  virtual void get_force_info(int, DL_actuator_type&,
			      DL_dyna*&, DL_point*, DL_vector*);

             DL_vtv();                       // constructor
	     ~DL_vtv();                      // destructor

  // for internal DL-use only:
  
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

  virtual void new_frame(void);
                     // prepare for a new frame
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
  virtual void get_error(DL_largevector*);
                     // calculate the constraint error vector
};

#endif
