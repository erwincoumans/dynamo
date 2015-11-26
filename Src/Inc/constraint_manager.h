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
// filename: constraint_manager.h
// description: system class for the inverse dynamics manager
// author: Bart Barenbrug   April 1996

#ifndef DL_CONSTRAINTMANAGERH
#define DL_CONSTRAINTMANAGERH

#include "list.h"
#include "largematrix.h"
#include "constraint.h"
#include "collision.h"

// the constraint manager administrates which constraint pairs
// influence each other, so it does not re-calculate zeros in
// the dCdR matrix all the time. Here is the class definition for
// those constrainmt-pair objects:

class DL_constraint_pair : public DL_ListElem {
  public:
    DL_constraint *cc;
    DL_constraint *cf;
    DL_constraint_pair(DL_constraint *c, DL_constraint *f):DL_ListElem()
    {  // constructor
      cc=c; cf=f;
    }
    ~DL_constraint_pair(){}; // destructor
}; // DL_constraint_pair
    

// *************************** //
// class DL_constraint_manager //
// *************************** //

class DL_constraint_manager {
  protected:
    DL_List	*c;          // The list of constraints that are managed
    DL_List cp;          // The list of constraint pairs that have
                         // a non-zero influence on each other
    int	dCdRToGo;        // number of frames to go before dCdR
                         // is recalculated
    int nr_cg;           // number of frames of conjug_grad solving with
                         // convergence
    boolean	c_changed;   // has the list of constraints changed since dCdR
                         // was calculated last
    DL_largematrix *dCdR;
    int	totdim;          // total dimension: sum of constraint dimensions
    DL_Scalar first_error;
    int nrcollisions;    // number of detected collisions
    DL_collision* *collisions; // since collisions are treated somewhat differently
    // references to them are listed here. I have to use an array to represent the list
    // since collisions are already present in the constraint list, and DL_ListElem-ents
    // can be in two lists at the same time.  B(
    int size_collisions;      // allocated size of the collisions-array
    boolean show_con_forces;
    
    void	new_frame(void);
                  // prepare all constraints for the new frame
    void	apply_first_estimates(void);
                  // for all constraints: calculate and apply first estimates
    void	reset_all(void);
                  // reset all constraints
    void	reset_undo_all(void);
                  // reset/undo all constraints
    void	calc_all_errors(DL_largevector *);
                  // calculates the large vector containing
                  // all coinstraint errors
    boolean	apply_all_restriction_changes(DL_largevector*);
                  // apply the restrictionchanges to all constraints
		  // returns if any constraints have deleted themselves in the process
    void	do_post_processing(void);
                  // give the constraints a change to do some post processing
    void    sort_constraints();
                  // based on the connectovity info in cp, sort the constraints so that
                  // the dCdR matrix has minimal bandwidth (or an approximation of that).
                  // also re-orders dCdR. This is useful for sparse LU decomposition.
    void	calc_dCdR_full(void);
    void	calc_dCdR_analytical(void);
    void	calc_dCdR_empirical(void);
    void    begin_test(void);
    void    end_test(void);
    void	iterate(DL_largevector*);
    void    redo_index_administration();
  public:
    /// control parameters etc. for external use:
    boolean     analytical; // analytical or empirical determination of
                            // dCdR
    int		MaxIter;    // the number of constraint correction iterations
    int		NrSkip;     // dCdR is recalculated every NrSkip+1 frames
    DL_Scalar   max_error;  // error thresh hold
    DL_Scalar   error;      // the current error magnitude
                            // (just for information)
    int         nriter;     // actual number of iteration steps taken in the
                            // last frame  (just for information)
    int         max_collisionloops;  // maximum number of secundary collision
                                     // detection/handling phases

    void    solve_using_lud();
    boolean solving_using_lud(){ return dCdR->get_solve_method()==lud_bcksub;};
    void    solve_using_cg();
    boolean solving_using_cg(){ return dCdR->get_solve_method()==conjug_grad;};
    void    solve_using_svd();
    boolean solving_using_svd(){ return dCdR->get_solve_method()==svd;};

    void    show_constraint_forces();
    void    hide_constraint_forces();
    boolean showing_constraint_forces() { return show_con_forces; };

    int     get_dof() { return totdim; };
    int     get_nr_constraints() { return c->length(); };

             DL_constraint_manager();           // constructor
	     ~DL_constraint_manager();          // destructor

    /// for internal DL use only:
    void	satisfy(void);       // do the constraint correction
    void	add(DL_constraint*); // add a constraint
    void	del(DL_constraint*); // delete a constraint

    void        add_collision(DL_collision*);
                // for DL_collision to be able to tell the constraint manager
		// that a collision has been detected (so secondary collision
		// detection can be called if necessary)

};

extern DL_constraint_manager* DL_constraints;

inline DL_constraint_manager::DL_constraint_manager() {
  if (DL_constraints && DL_constraints!=this)
    DL_dsystem->get_companion()->Msg("Error: there should only be one constraint manager!!\n");
  if (!DL_constraints) DL_constraints=this;
  MaxIter=10;
  NrSkip=dCdRToGo=totdim=0;
  error=first_error=0.0; max_error=0.1;
  nriter=nr_cg=0;
  max_collisionloops=1; // no secundary collision detection by default
  c_changed=FALSE;
  analytical=TRUE;
  c=new DL_List;
  dCdR=new DL_largematrix(0,0);
  nrcollisions=size_collisions=0;
  show_con_forces=FALSE;
}

inline DL_constraint_manager::~DL_constraint_manager(void) {
  cp.delete_all();
  if (size_collisions>0) delete[] collisions;
  delete dCdR;
  delete c;
}

inline void DL_constraint_manager::redo_index_administration(){
  totdim=0;
  DL_constraint *celem=(DL_constraint *)c->getfirst();
  while (celem) {
    celem->index=totdim;
    totdim+=celem->dim;
    celem=(DL_constraint*)c->getnext(celem);
  }
}

inline void DL_constraint_manager::new_frame(void) {
// note that constraints can delete themselves in this phase, so
// some extra care is required in traversing c:
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  DL_constraint *next_c;
  while (constr) {
    next_c=(DL_constraint*)c->getnext(constr);
    constr->new_frame();
    constr=next_c;
  }
}

inline void DL_constraint_manager::apply_first_estimates(void) {
// note that constraints can delete themselves in this phase, so
// some extra care is required in traversing c:
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  DL_constraint *next_c;
  while (constr) {
    next_c=(DL_constraint*)c->getnext(constr);
    constr->first_estimate();
    constr=next_c;
  }
}

inline void DL_constraint_manager::reset_all(void) {
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  while (constr) {
    constr->reset();
    constr=(DL_constraint*)c->getnext(constr);
  }
}

inline void DL_constraint_manager::reset_undo_all(void) {
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  while (constr) {
    constr->reset_undo();
    constr=(DL_constraint*)c->getnext(constr);
  }
}

inline void DL_constraint_manager::do_post_processing(void) {
// note that constraints can delete themselves in this phase, so
// some extra care is required in traversing c:
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  DL_constraint *next_c;
  while (constr) {
    next_c=(DL_constraint*)c->getnext(constr);
    constr->post_processing();
    constr=next_c;
  }
}

inline void DL_constraint_manager::begin_test() {
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  while (constr) {
    constr->begin_test();
    constr=(DL_constraint*)c->getnext(constr);
  }
}

inline void DL_constraint_manager::end_test() {
  DL_constraint *constr=(DL_constraint *)c->getfirst();
  while (constr) {
    constr->end_test();
    constr=(DL_constraint*)c->getnext(constr);
  }
}

#endif
