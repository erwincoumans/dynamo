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
// filename: pts.h
// description: point-to-surface constraint class
// author: Bart Barenbrug   June '96
//

#ifndef DL_PTSH
#define DL_PTSH

#include "constraint.h"
#include "constraint_manager.h"
#include "surface.h"

// ************ //
// class DL_pts //
// ************ //

class DL_pts : public DL_constraint {
 protected:
  DL_surface *surf;   // the constraint connects surf(s,t) to p in g.
  DL_Scalar s,olds,ssave; // either g or c->get_geo() is a dyna.
  DL_Scalar t,oldt,tsave;
  DL_Scalar sf,tf; // surface parameters for point to which forces are applied
  DL_geo   *g;
  DL_point p;
  DL_point pgw;   // p in wc if !g_is_dyna
  DL_vector dpgw; // velocity of pg (in wc; used only if !g_is_dyna)
  boolean g_is_dyna, sg_is_dyna;

  DL_vector x,y,n; // the derivatives for the current frame and the normal that
                   // make the local local coordinate system (expressed in
	  	   // world coordinates).
  DL_Scalar xxinv;     // 1.0/x->inprod(x)
  DL_Scalar yyinv;     // 1.0/y->inprod(y)
  DL_largematrix dcdp; // the matrix dc/dp of which the one row consist of n
  DL_largematrix dfdR; // the matrix df/dR of which the one column consist of n
  DL_point sstl; // the surface position for this frame (in local coordinates);
                 // only valid if sg_is_dyna;  
  DL_point sstf; // surface(sf,tf) (in local coordinates)
  boolean st_inbounds;
  
  void getforce(DL_largevector*, DL_vector*);
                     // given the given restrictions return the corresponding force
public:
  /// for external use:
  DL_Scalar maxforce; // maximum force (<=0: no maxforce checking)

  void init(DL_geo*, DL_point*, DL_surface*); // initialise the constraint
  void reactionforce(DL_vector*);
                     // returns a copy of the reaction force (in world
		     // coordinates)
  DL_geo* get_geo(void){ return g; };
  void get_point(DL_point *pnt){ pnt->assign(&p); };
  void set_point(DL_point *pnt){ p.assign(pnt); };
  DL_surface *get_surface(void){ return surf; };
  DL_Scalar get_s(void); // get current first surface parameter
  DL_Scalar get_t(void); // get current second surface parameter

  virtual void get_fd_info(int&,int&);
  virtual void get_force_info(int, DL_actuator_type&,
			      DL_dyna*&, DL_point*, DL_vector*);

             DL_pts();                    // constructor
	     ~DL_pts(){};                 // destructor

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
