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
// filename: linehinge.h
// description: linehinge constraint class
// author: Bart Barenbrug   July '96
//

#ifndef DL_LINEHINGEH
#define DL_LINEHINGEH

#include "constraint_manager.h"

// ****************** //
// class DL_linehinge //
// ****************** //

class DL_linehinge : public DL_constraint {
protected:
  DL_dyna  *d;         // the constraint connects points pd0 and pd1 in d to
  DL_point pd0,pd1;  // point pg0 and pg1 in g respectively
  DL_geo   *g;
  DL_point pg0,pg1;
  DL_point pg0w,pg1w; // pg? in wc (only calculated once if !g_is_dyna)
  DL_vector dpg0w,dpg1w; // velocity of pg? (in wc; used only if !g_is_dyna) 
  boolean g_is_dyna;
  DL_vector l,x,y; // local coordinate system (in wc)
  DL_largematrix dcdp0, dcdp1; // dc/dp matrices
  DL_largematrix df0dR, df1dR; // df/dR matrices

  void getforce0(DL_largevector*,DL_vector*);  // transform restriction vector to forces:
  void getforce1(DL_largevector*,DL_vector*);  // 0: in  p?0; 1: in p?1
  void getforce2(DL_largevector*,DL_vector*);  // 2: along l;

public:
  /// externally accessible:
  DL_Scalar maxforce; // maximum force magnitude:
                  //   (<=0: no maxforce checking)
		  // the maximum is taken in terms of the length of the
		  // 5D-restriction value
  void init(DL_dyna*, DL_point*, DL_point*, DL_geo*, DL_point*, DL_point*);
       // initialise the constraint

  DL_dyna* get_dyna(void){ return d; };
  void get_dyna_point0(DL_point *dp){ dp->assign(&pd0); };
  void set_dyna_point0(DL_point *dp){ pd0.assign(dp); };
  void get_dyna_point1(DL_point *dp){ dp->assign(&pd1); };
  void set_dyna_point1(DL_point *dp){ pd1.assign(dp); };
  DL_geo* get_geo(void){ return g; };
  void get_geo_point0(DL_point *gp){ gp->assign(&pg0); };
  void set_geo_point0(DL_point *gp){ pg0.assign(gp); };
  void get_geo_point1(DL_point *gp){ gp->assign(&pg1); };
  void set_geo_point1(DL_point *gp){ pg1.assign(gp); };

  void reactionforce0(DL_vector*);  // the reaction forces:
  void reactionforce1(DL_vector*);  // 0: in  p?0; 1: in p?1
  void reactionforce2(DL_vector*);  // 2: along l;

  DL_Scalar reactionforce(void);  // length of current restriction vector;

  virtual void get_fd_info(int&,int&);
  virtual void get_force_info(int, DL_actuator_type&,
			      DL_dyna*&, DL_point*, DL_vector*);

             DL_linehinge();                    // constructor
	     ~DL_linehinge();                   // destructor

  /// for DL-internal use only:
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
  virtual void post_processing(void);
                     // post processing

};

#endif
