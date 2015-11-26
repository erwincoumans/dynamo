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
// filename     : constraint.cpp
// description	: non-inline methods of class DL_constraint
// author	: Bart Barenbrug     April 1996
//

//#define DEBUG

#include "constraint.h"
#include "constraint_manager.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_constraint::DL_constraint():DL_ListElem(),DL_force_drawable() {
  index=0;
  stiffness=1.0;
  dim=0;
  F=new DL_largevector(dim);
  oldF=new DL_largevector(dim);
  Fsave=new DL_largevector(dim);
  active=initialised=testing=FALSE;
  veloterms=FALSE;
  veloterms_free=TRUE;
  nr_osc=4;
  max_osc=8;
}

DL_constraint::~DL_constraint(void) {
  deactivate();
  delete Fsave;
  delete oldF;
  delete F;
}

void DL_constraint::init(void) {
  initialised=TRUE;
  activate();
  nr_osc=(rand()%max_osc);
}

void DL_constraint::rotatebase(DL_vector *a, DL_vector *b, DL_vector *x, DL_vector *y){
// PRE: a=A && b=B && x=X && y=Y
// POST: a=RA && x=RX && y=RY where R is a rotation satisfying RA=B
// (the rotation around crossprod(A,B) is used)

// We use local coordinate system <a, a.cross(b), a.cross(a.crossprod(b))>
// of which the rotated version is:
// <b, a.crossprod(b), b.crossprod(a.crossprod(b))>
// Vectors x and y are first expressed in this local coordinate system
// enabling the use of the local coordinates in the rotated version of the
// coordinate system to arrive at the rotated versions of x and y.

  DL_vector anorm,bnorm,ab,aab,bab,z;
  DL_Scalar za, zab, zaab;
  
  anorm.assign(a); anorm.normalize();
  bnorm.assign(b); bnorm.normalize();
  anorm.crossprod(&bnorm,&ab);
  ab.normalize();
  if (ab.norm()<0.9) return; // no rotation required
  anorm.crossprod(&ab,&aab);
  bnorm.crossprod(&ab,&bab);

  // now rotate x:
  za=x->inprod(&anorm);
  zab=x->inprod(&ab);
  zaab=x->inprod(&aab);
  bnorm.times(za,x);
  ab.times(zab,&z);
  x->plusis(&z);
  bab.times(zaab,&z);
  x->plusis(&z);

  // now rotate y:
  za=y->inprod(&anorm);
  zab=y->inprod(&ab);
  zaab=y->inprod(&aab);
  bnorm.times(za,y);
  ab.times(zab,&z);
  y->plusis(&z);
  bab.times(zaab,&z);
  y->plusis(&z);

  // and finaly rotate a:
  a->assign(b);
}


void DL_constraint::soft(void) {
  veloterms=TRUE;
  veloterms_free=FALSE;
}

void DL_constraint::hard(void) {
  veloterms=FALSE;
  veloterms_free=FALSE;
}

void DL_constraint::auto_softhard() {
  veloterms_free=TRUE;
  nr_osc=(index%max_osc);
}

void DL_constraint::auto_softhard(int i) {
  max_osc=i;
  veloterms_free=TRUE;
  nr_osc=(index%max_osc);
}

void DL_constraint::activate(void) {
  if (active) return;
  if (!initialised) {
    DL_dsystem->get_companion()->Msg("Cannot activate an uninitialised constraint!\n");
    return;
  }
  if (DL_constraints) {
     DL_constraints->add(this);
     reset();
     active=TRUE;
  }
  else
     DL_dsystem->get_companion()->Msg("Can't activate constraint because there is no constraint manager!\n");
  
}

void DL_constraint::deactivate(void) {
  if (active) {
    if (DL_constraints) DL_constraints->del(this);
    else DL_dsystem->get_companion()->Msg("constraint::deactivate(): no constraint manager!\n");
    reset();
  }
  active=FALSE;
}

void DL_constraint::begin_test(void) {
  Fsave->assign(F);
  testing=TRUE;
}

void DL_constraint::end_test(void) {
  F->assign(Fsave);
  testing=FALSE;
}

boolean DL_constraint::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
#ifdef DEBUG
  DL_dsystem->get_companion()->Msg("constraint::dCdRsub called!!\n");
#endif
  return FALSE;
}

boolean DL_constraint::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
#ifdef DEBUG
  DL_dsystem->get_companion()->Msg("constraint::dCdFq called!!\n");
#endif
  return FALSE;
}

boolean DL_constraint::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
#ifdef DEBUG
  DL_dsystem->get_companion()->Msg("constraint::dCdF called!!\n");
#endif
  return FALSE;
}

boolean DL_constraint::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
#ifdef DEBUG
  DL_dsystem->get_companion()->Msg("constraint::dCdM called!!\n");
#endif
  return FALSE;
}

boolean DL_constraint::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdi has dimensions dim x 3
#ifdef DEBUG
  DL_dsystem->get_companion()->Msg("constraint::dCdI called!!\n");
#endif
  return FALSE;
}

void DL_constraint::new_frame(void) {
  detect_osc();
}

void DL_constraint::first_estimate(void) {
  static DL_largevector dF(dim);
  if (oldF->norm()==0.0) {
    // constant extrapolation:
    oldF->assign(F);
  }
  else {
    // linear extrapolation:
    F->minus(oldF,&dF);
    oldF->assign(F);
    dF.timesis(stiffness);
    F->plusis(&dF);
  }
//   if (check_restrictions()) // removed since first estimate is too much of a
                               // guess and can be way too high (esp. after
			       // initialization)
   apply_restrictions(F);
}

boolean DL_constraint::check_restrictions() {
  return TRUE;
}

void DL_constraint::apply_restrictions(DL_largevector* lv) {
// nothing to do for empty constraint
#ifdef DEBUG
  DL_dsystem->get_companion()->Msg("constraint::apply_restrictions called!!\n");
#endif
}

void DL_constraint::test_restriction_changes(DL_largevector* lv) {
// nothing to do for an empty constraint
}

void DL_constraint::apply_restriction_changes(DL_largevector* lv) {
  F->plusis(lv);
  apply_restrictions(lv);
}

void DL_constraint::get_error(DL_largevector* lv) {
// nothing to do for empty constraint
}

void DL_constraint::detect_osc(void) {
// set veloterms if there are any oscilations detected.
// to be called before constraint::first_estimate() in the
//   first_estimate method of any constraint specialisation
//   with positional terms

  if (veloterms_free) {
    veloterms=(nr_osc==0);
    nr_osc++; if (nr_osc>=max_osc) nr_osc=0;
  }
  
// maybe add:
// if (veloterms) { /* raise an event */ }
}

void DL_constraint::post_processing(void) {
  // nothing to do for empty constraint
}

void DL_constraint::reset_undo(void){
  F->neg(F);
  apply_restrictions(F);
  reset();
}
#undef DEBUG
