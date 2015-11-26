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
// filename	: supvec.cpp
// description	: non-inline methods of class DL_supvec
// author	: Bart Barenbrug   May '96
//

#include "dyna.h"
#include "supvec.h"

// *************** //
// member fuctions //
// *************** //

DL_supvec::DL_supvec(){
  z.init(0,0,0);
  v.init(0,0,0);
  q.init(0,0,0,1);
  A.assign(1,0,0, 0,1,0, 0,0,1);
  w.init(0,0,0);
}

DL_supvec::DL_supvec(DL_point* nz, DL_vector* nv,
                            DL_vector4* nq, DL_vector* nw, DL_matrix* nA) {
  z.assign(nz); v.assign(nv); q.assign(nq); w.assign(nw); A.assign(nA);
}

void DL_supvec::init(DL_point* nz, DL_vector* nv,
                            DL_vector4* nq, DL_vector* nw, DL_matrix *nA) {
  z.assign(nz); v.assign(nv); q.assign(nq); w.assign(nw); A.assign(nA);
}

void DL_supvec::assign(DL_supvec* sv) {
  if (sv) {
    z.assign(&(sv->z)); v.assign(&(sv->v));
    q.assign(&(sv->q)); w.assign(&(sv->w));
    A.assign(&(sv->A));
  }
}

// the following methods only operate on q and w since they are used by
// the motion integrators of dyna

void DL_supvec::plusis(DL_supvec* sv) { 
  if (sv) {
    q.plusis(&(sv->q));
    w.plusis(&(sv->w));
  }
}

void DL_supvec::timesis(DL_Scalar f) { 
  q.timesis(f);
  w.timesis(f);
}

void DL_supvec::plus(DL_supvec* sv,DL_supvec *nv) { 
  if (sv) {
    q.plus(&(sv->q), &(nv->q));
    w.plus(&(sv->w), &(nv->w));
  }
}

void DL_supvec::times(DL_Scalar t, DL_supvec *nv) {
  q.times(t, &(nv->q));
  w.times(t, &(nv->w));
}

void DL_supvec::A2q() {
  q.from_matrix(&A);
}

void DL_supvec::q2A(DL_dyna *d) {
// q.normalize();
 q.to_matrix(&A);

 // for 1D-objects the angular velocity should not have any components in the
 // direction of the object's axis. Since we changed the orientation of that
 // axis, we might have to adjust the angular velocity as well:
  if (d->get_oneD()!=0) {
    DL_vector *l;
    if (d->get_oneD()==1) l=&(A.c0);
    if (d->get_oneD()==2) l=&(A.c1);
    if (d->get_oneD()==3) l=&(A.c2);
    DL_vector wp;
    l->times(l->inprod(&w),&wp);
    // wp is projection of w onto the axis of the object: this
    // is what we want to get rid of:
    w.minus(&wp,&w);
  }
}
