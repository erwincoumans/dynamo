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
// filename	: geo.h
// description	: geometry companion objects
// author	: Bart Barenbrug   Januari '97
//

#ifndef DL_GEOH
#define DL_GEOH

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include "pointvector.h"
#include "boolean.h"
#include "supvec.h"
#include "list.h"

// ************ //
// class DL_geo //
// ************ //

class DL_geo : public DL_ListElem {
  protected:
    void *companion; // the companion from the user system
    // motion state:
    DL_supvec mstate;  // the motion state at time t: q is not used (A is)
    DL_supvec nextmstate;  // the motion state at t+h: q is not used (A is)

    // material properties:
    DL_Scalar elasticity;
    
  public:
    void* get_companion(){ return companion; };
    virtual void set_position(DL_point*);     // set the geo's position
    virtual DL_point* get_position(void);     // get the geo's position
    virtual void set_velocity(DL_vector*);    // set the geo's velocity
    virtual DL_vector* get_velocity(void);    // get the geo's velocity
    virtual void set_orientation(DL_matrix*); // set the geo's orientation
    virtual DL_matrix* get_orientation(void); // get the geo's orientation
    virtual void set_angvelocity(DL_vector*); // set the geo's angular velocity
    virtual DL_vector* get_angvelocity(void); // get the geo's angular velocity

    virtual void move(DL_point*,DL_matrix*);  // give the geo the new position
                                              // and orientation given by the
			                      // parameters, and update the
				              // (angular) velocity accordingly

    virtual void to_world(DL_point*,DL_point*);
    virtual void to_world(DL_vector*,DL_vector*);
      // converts the coordinates of the first point/vector (given in local
      // coordinates) to world coordinates based on mstate
    virtual void to_local(DL_point*,DL_geo*,DL_point*);
    virtual void to_local(DL_vector*,DL_geo*,DL_vector*);
      // convert the coordinates of the first point/vector (as expressed in
      // the coordinate system of the DL_geo, or in world coordinates if that
      // pointer is zero) to local coordinates (based on mstate)
    virtual void get_velocity(DL_point*,DL_vector*);
    virtual void get_velocity(DL_vector*,DL_vector*);
      // calculates the velocity (in wc) of the given (in lc) point/vector
      // based on mstate

    void assign(DL_geo*,void*); // assignment (needs reference to new
                                // companion)

    virtual boolean is_dyna(void) { return FALSE; }

    void  set_elasticity(DL_Scalar el){ elasticity=el;};
    DL_Scalar get_elasticity(void){ return elasticity;};
  
    DL_geo(void*); //constructor
    ~DL_geo(){};

    /// For internal use (by the constraints)::
  
    virtual void set_next_position(DL_point*);     // set the geo's next position
    virtual DL_point* get_next_position(void);     // get the geo's next position
    virtual void set_next_velocity(DL_vector*);    // set the geo's next velocity
    virtual DL_vector* get_next_velocity(void);    // get the geo's next velocity
    virtual void set_next_orientation(DL_matrix*); // set the geo's next orientation
    virtual DL_matrix* get_next_orientation(void); // get the geo's next orientation
    virtual void set_next_angvelocity(DL_vector*); // set the geo's next angular velocity
    virtual DL_vector* get_next_angvelocity(void); // get the geo's next angular velocity

    virtual void new_toworld(DL_point*, DL_point*);
    virtual void new_toworld(DL_vector*,DL_vector*);
      // converts the coordinates of the first point/vector (given in local
      // coordinates) to world coordinates based on nextmstate
    virtual void new_tolocal(DL_point*,DL_geo*,DL_point*);
    virtual void new_tolocal(DL_vector*,DL_geo*,DL_vector*);
      // convert the coordinates of the first point/vector (as expressed in
      // the coordinate system of the DL_geo, or in world coordinates if that
      // pointer is zero) to local coordinates (based on nextmstate)
    virtual void get_newvelocity(DL_point*, DL_vector*);
    virtual void get_newvelocity(DL_vector*,DL_vector*);
      // calculates the velocity (in wc) of the given (in lc) point/vector
      // based on nextmstate
};

inline DL_geo::DL_geo(void* comp) : DL_ListElem() {
  companion=comp;
  elasticity=1.0;
}

inline void DL_geo::assign(DL_geo *g, void *newcomp){
    companion=newcomp;
    mstate.assign(&(g->mstate));
}

inline void DL_geo::set_position(DL_point *p) {
  mstate.z.assign(p);
}

inline DL_point* DL_geo::get_position(void) {
  return &(mstate.z);
}

inline void DL_geo::set_velocity(DL_vector *v) {
  mstate.v.assign(v);
}

inline DL_vector* DL_geo::get_velocity(void) {
  return &(mstate.v);
}

inline void DL_geo::set_orientation(DL_matrix *m) {
  mstate.A.assign(m);
}

inline DL_matrix* DL_geo::get_orientation(void) {
  return &(mstate.A);
}

inline void DL_geo::set_angvelocity(DL_vector *w) {
  mstate.w.assign(w);
}

inline DL_vector* DL_geo::get_angvelocity(void) {
  return &(mstate.w);
}

inline void DL_geo::set_next_position(DL_point *p) {
  nextmstate.z.assign(p);
}

inline DL_point* DL_geo::get_next_position(void) {
  return &(nextmstate.z);
}

inline void DL_geo::set_next_velocity(DL_vector *v) {
  nextmstate.v.assign(v);
}

inline DL_vector* DL_geo::get_next_velocity(void) {
  return &(nextmstate.v);
}

inline void DL_geo::set_next_orientation(DL_matrix *m) {
  nextmstate.A.assign(m);
}

inline DL_matrix* DL_geo::get_next_orientation(void) {
  return &(nextmstate.A);
}

inline void DL_geo::set_next_angvelocity(DL_vector *w) {
  nextmstate.w.assign(w);
}

inline DL_vector* DL_geo::get_next_angvelocity(void) {
  return &(nextmstate.w);
}

inline void DL_geo::to_world(DL_point *pl, DL_point *pw) {
  DL_vector vtmp;
  mstate.A.times(pl,pw);
  mstate.z.tovector(&vtmp);
  pw->plusis(&vtmp);
}

inline void DL_geo::to_world(DL_vector *vl, DL_vector *vw) {
  mstate.A.times(vl,vw);
}

inline void DL_geo::to_local(DL_point *p, DL_geo *g, DL_point *pl) {
  if (g==this) {
    pl->assign(p);
    return;
  }
  if (g==NULL) {
    DL_vector vtmp,vl;
    DL_matrix Ainv;
    p->minus(&(mstate.z),&vtmp);
    mstate.A.invert(&Ainv);
    Ainv.times(&vtmp,&vl);
    vl.topoint(pl);
    return;
  }
  DL_point pw;
  g->to_world(p,&pw);
  to_local(&pw,NULL,pl);
}

inline void DL_geo::to_local(DL_vector *v, DL_geo *g, DL_vector *vl) {
  if (g==this) {
    vl->assign(v);
    return;
  }
  if (g==NULL) {
    DL_matrix Ainv;
    mstate.A.invert(&Ainv);
    Ainv.times(v,vl);
    return;
  }
  DL_vector vw;
  g->to_world(v,&vw);
  to_local(&vw,NULL,vl);
}

inline void DL_geo::get_velocity(DL_point *p, DL_vector *v){
  // p in local coordinates; v in world coordinates
  DL_point pw;
  DL_vector vw;
  mstate.A.times(p,&pw);
  pw.tovector(&vw);
  vw.crossprod(&(mstate.w),v);
  v->plusis(&(mstate.v));  
}

inline void DL_geo::get_velocity(DL_vector *v, DL_vector *vr){
  DL_vector vw;
  mstate.A.times(v,&vw);
  vw.crossprod(&(mstate.w),vr);
}

inline void DL_geo::new_toworld(DL_point *p, DL_point *pr) {
  DL_vector vtmp;
  nextmstate.A.times(p,pr);
  nextmstate.z.tovector(&vtmp);
  pr->plusis(&vtmp);
}

inline void DL_geo::new_toworld(DL_vector *v, DL_vector *vr) {
  nextmstate.A.times(v,vr);
}

inline void DL_geo::new_tolocal(DL_point *p, DL_geo *g, DL_point *pl) {
  if (g==this) {
    pl->assign(p);
    return;
  }
  if (g==NULL) {
    DL_vector vtmp,vl;
    p->minus(&(nextmstate.z),&vtmp);
    nextmstate.A.transposetimes(&vtmp,&vl); // A's inverse is its transpose
    vl.topoint(pl);
    return;
  }
  DL_point pw;
  g->new_toworld(p,&pw);
  new_tolocal(&pw,NULL,pl);
}

inline void DL_geo::new_tolocal(DL_vector *v, DL_geo *g, DL_vector *vl) {
  if (g==this) {
    vl->assign(v);
    return;
  }
  if (g==NULL) {
    nextmstate.A.transposetimes(v,vl); // A's inverse is its transpose
    return;
  }
  DL_vector vw;
  g->new_toworld(v,&vw);
  new_tolocal(&vw,NULL,vl);
}

inline void DL_geo::get_newvelocity(DL_point *p, DL_vector *v) {
// p in local coordinates; v in world coordinates
  DL_point pw;
  DL_vector vw;
  nextmstate.A.times(p,&pw);
  pw.tovector(&vw);
  vw.crossprod(&(nextmstate.w),v);
  v->plusis(&(nextmstate.v));
}

inline void DL_geo::get_newvelocity(DL_vector *v, DL_vector *vr) {
// v in local coordinates; vr in world coordinates
  DL_vector vw;
  nextmstate.A.times(v,&vw);
  vw.crossprod(&(nextmstate.w),vr);
}

#endif
