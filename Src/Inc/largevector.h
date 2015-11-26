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
// filename	: largevector.h
// description	: arbitrary-D vectors (dimension is known in advance)
// author	: Bart Barenbrug   April '96
//

#ifndef DL_LARGEVECTORH
#define DL_LARGEVECTORH

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include "pointvector.h"
#include "boolean.h"
#include "scalar.h"

// ******************** //
// class DL_largevector //
// ******************** //

class DL_largevector {
  friend class DL_largematrix;
  protected:
    int        dim;  // the dimension of the vector
    int        vsize; // size of v-array; dim<=vsize:
    DL_Scalar* v; // the coordinate array

  public:
    int        get_dim() { return dim; };
    void       init(DL_Scalar);
    void       init(DL_Scalar,DL_Scalar);
    void       init(DL_Scalar,DL_Scalar,DL_Scalar);
    void       init(DL_Scalar,DL_Scalar,DL_Scalar,DL_Scalar);
    void       init(DL_Scalar,DL_Scalar,DL_Scalar,DL_Scalar,DL_Scalar);
    void       init(DL_Scalar,DL_Scalar,DL_Scalar,DL_Scalar,DL_Scalar,DL_Scalar);
    void       resize(int);
    void       makezero(void);
    
    DL_Scalar  get(int);
    void       set(int,DL_Scalar);
    void       setsubvector(int, DL_largevector*);
    void       getsubvector(int, DL_largevector*);

    void       assign(DL_largevector*);
    void       assign(DL_vector*);
    DL_Scalar  norm();
    void       normalize();
    DL_Scalar  inprod(DL_largevector*);
    void       plusis(DL_largevector*);
    void       minusis(DL_largevector*);
    void       timesis(DL_Scalar);
    boolean    equal(DL_largevector*);

    void       neg(DL_largevector*);
    void       times(DL_Scalar,DL_largevector*);       
    void       plus(DL_largevector*,DL_largevector*);      
    void       minus(DL_largevector*,DL_largevector*);
    
          DL_largevector(int dimension=0);  // constructor
          DL_largevector(DL_largevector*);  // copy constructor
	  ~DL_largevector();                // destructor

    // for debugging:
    void       show(void); // print the vector using Msg()
};

inline DL_largevector::DL_largevector(int dimension) {
  vsize=dim=dimension;
  if (vsize<6) vsize=6;
  v=new DL_Scalar[vsize];
}

inline void DL_largevector::resize(int newdim) {
// This loses all info stored in the vector!!!
  if (dim==newdim) return;
  dim=newdim;
  if ((dim>vsize) || ((dim>6) && (2*dim<vsize))) {
    if (v) delete[] v;
    vsize=dim;
    v=new DL_Scalar[vsize];
  }
}

inline void DL_largevector::assign(DL_largevector* lv) {
// PRE: lv
  resize(lv->dim);
  for (register int i=0;i<dim;i++) v[i]=lv->v[i];
}

inline DL_largevector::DL_largevector(DL_largevector* lv) {
  v=NULL;
  assign(lv);
}

inline DL_largevector::~DL_largevector() {
  delete[] v;
}

inline void DL_largevector::init(DL_Scalar c0) {
  if (dim>=1) v[0]=c0;
}

inline void DL_largevector::init(DL_Scalar c0, DL_Scalar c1) {
  if (dim>=2) {
     v[0]=c0; v[1]=c1;
  }
}

inline void DL_largevector::init(DL_Scalar c0, DL_Scalar c1, DL_Scalar c2) {
  if (dim>=3) {
     v[0]=c0; v[1]=c1; v[2]=c2;
  }
}

inline void DL_largevector::init(DL_Scalar c0, DL_Scalar c1, DL_Scalar c2, DL_Scalar c3) {
  if (dim>=4) {
     v[0]=c0; v[1]=c1; v[2]=c2; v[3]=c3;
  }
}

inline void DL_largevector::init(DL_Scalar c0, DL_Scalar c1, DL_Scalar c2, DL_Scalar c3, DL_Scalar c4) {
  if (dim>=5) {
     v[0]=c0; v[1]=c1; v[2]=c2; v[3]=c3; v[4]=c4;
  }
}

inline void DL_largevector::init(DL_Scalar c0, DL_Scalar c1, DL_Scalar c2, DL_Scalar c3, DL_Scalar c4, DL_Scalar c5) {
  if (dim>=6) {
     v[0]=c0; v[1]=c1; v[2]=c2; v[3]=c3; v[4]=c4; v[5]=c5;
  }
}

inline void DL_largevector::assign(DL_vector* vec) {
// PRE: dim>=3 && vec
  v[0]=vec->x; v[1]=vec->y; v[2]=vec->z;
}

inline void DL_largevector::makezero(void) {
  for(register int i=0;i<dim;i++) v[i]=0;
}

inline DL_Scalar DL_largevector::norm() {
  register DL_Scalar l=0.0;
  register int i;
  for (i=0; i<dim; i++) l+=v[i]*v[i];
  return ::sqrt(l);
}

inline void DL_largevector::normalize() {
  register DL_Scalar l=norm();
  if (l!=0) {
    l=1.0/l;
    for (register int i=0; i<dim; i++) v[i]*=l;
  }
}

inline DL_Scalar DL_largevector::inprod(DL_largevector* lv) {
// PRE: (lv->dim>=dim)
  double inp=0.0;
  if (lv) {
    for (register int i=0; i<dim ;i++) inp+=v[i]*lv->v[i];
  }
  return (DL_Scalar)inp;
}

inline void DL_largevector::plusis(DL_largevector* lv) {
// PRE: lv && lv->dim>=dim
  for (register int i=0;i<dim;i++) v[i]+=lv->v[i];
}

inline void DL_largevector::minusis(DL_largevector* lv) {
// PRE: lv && lv->dim>=dim
  for (register int i=0;i<dim;i++) v[i]-=lv->v[i];
}

inline void DL_largevector::timesis(DL_Scalar f) {
   for (register int i=0; i<dim; i++) v[i]*=f;
}

inline boolean DL_largevector::equal(DL_largevector* lv) {
  if (!lv) return FALSE;
  if (dim!=lv->dim) return FALSE;
  register int i=0;
  while ((i<dim) && (v[i]==lv->v[i])) i++;
  return (i==dim);
}

inline DL_Scalar DL_largevector::get(int r) {
  if ((0<=r) && (r<dim)) return v[r];
  else return 0;
}

inline void DL_largevector::set(int r,DL_Scalar f) {
  if ((0<=r) && (r<dim)) v[r]=f;
}

inline void DL_largevector::plus(DL_largevector* lv, DL_largevector *nlv) {
// PRE: lv && nlv && (lv->dim>=dim) && (nlv->dim>=dim)
  for(register int i=0;i<dim;i++) nlv->v[i]=v[i]+lv->v[i];
}

inline void DL_largevector::minus(DL_largevector* lv, DL_largevector *nlv) {
// PRE: lv && nlv && (lv->dim>=dim) && (nlv->dim>=dim)
  for(register int i=0;i<dim;i++) nlv->v[i]=v[i]-lv->v[i];
}

inline void DL_largevector::neg(DL_largevector *nlv) {
// PRE: nlv && (nlv->dim>=dim)
  for(register int i=0;i<dim;i++) nlv->v[i]=-v[i];
}

inline void DL_largevector::times(DL_Scalar t, DL_largevector *nlv) {
// PRE: nlv && (nlv->dim>=dim)
  for(register int i=0;i<dim;i++) nlv->v[i]=v[i]*t;
}

inline void DL_largevector::setsubvector(int idx, DL_largevector* lv) {
//PRE: lv && (idx+lv->dim<=dim)
  register int idxi=idx;  // INV: idxi==idx+i
  for (register int i=0; i<lv->dim; i++) { v[idxi]=lv->v[i]; idxi++; }
}

inline void DL_largevector::getsubvector(int idx, DL_largevector* lv) {
//PRE: lv && (idx+lv->dim<=dim)
  register int idxi=idx;  // INV: idxi==idx+i
  for (register int i=0; i<lv->dim; i++) { lv->v[i]=v[idxi]; idxi++; }
}

#endif
