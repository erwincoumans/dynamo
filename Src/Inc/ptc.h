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
// filename: ptc.h
// description: point-to-curve constraint class
// author: Bart Barenbrug   June '96
//

#ifndef DL_PTCH
#define DL_PTCH

#include "constraint.h"
#include "constraint_manager.h"
#include "curve.h"

// ************ //
// class DL_ptc //
// ************ //

class DL_ptc : public DL_constraint {
 protected:
  DL_curve  *c;           // the constraint connects c(s) to p in g.
  DL_Scalar s,olds,ssave; // either g or c->get_geo() is a dyna.
  DL_Scalar sf;           // curve parameter for point to which forces are applied
  DL_geo   *g;
  DL_point p;
  DL_point pgw;    // p in wc if !g_isdyna
  DL_vector dpgw;  // velocity of pg (in wc; used only if !g_is_dyna)
  boolean g_is_dyna, cg_is_dyna;

  DL_vector d,x,y; // the derivative for the current frame and the two vectors that
                   // make the local local coordinate system with d (expressed in
		   // world coordinates).
  DL_Scalar ddinv;     // 1.0/d->inprod(d)
  DL_largematrix dcdp; // the matrix dc/dp of which the two rows consist of
                       // x and y
  DL_largematrix dfdR; // the matrix df/dR of which the two columns consist of
                       // x and y
  DL_point csl;   // the curve position for this frame (in local coordinates);
                   // only valid if cg_is_dyna;
  DL_point csf;    // curve(sf) (local coordinates)
  boolean s_inbounds;
  
  void getforce(DL_largevector*, DL_vector*);
                     // given the given restrictions return the corresponding force
public:
  /// for external use:
  DL_Scalar maxforce; // maximum force (<=0: no maxforce checking)

  void init(DL_geo*, DL_point*, DL_curve*); // initialise the constraint
  void reactionforce(DL_vector*);
                     // returns a copy of the reaction force (in world
		     // coordinates)
  DL_geo* get_geo(void){ return g; };
  void get_point(DL_point *pnt){ pnt->assign(&p); };
  void set_point(DL_point *pnt){ p.assign(pnt); };
  DL_curve *get_curve(void){ return c; };
  DL_Scalar get_s(void); // get current curve parameter
		     
  virtual void get_fd_info(int&,int&);
  virtual void get_force_info(int, DL_actuator_type&,
			      DL_dyna*&, DL_point*, DL_vector*);

             DL_ptc();                    // constructor
	     ~DL_ptc(){};                 // destructor

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
  

  virtual void new_frame(void);
                     // prepare for the new frame
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
