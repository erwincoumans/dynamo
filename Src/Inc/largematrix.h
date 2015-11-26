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
// filename	: largematrix.h
// description	: arbitrary-D matrices (dimensions known in advance)
//                These are geared toward regular use as full matrices
//                and use for solving Ax=b where the square matrix A is
//                possibly sparse.
// author	: Bart Barenbrug   Sept 1997
//

#ifndef LARGEMATRIXH
#define LARGEMATRIXH

#include "boolean.h"
#include "matrix.h"
#include "largevector.h" 
#include "minmax.h"
#include "dyna_system.h"
enum solve_method {lud_bcksub, conjug_grad, svd};

// ******************** //
// class DL_largematrix //
// ******************** //

class DL_largematrix {
  protected:

    enum representation {full, riss, lud, ludb, svdcmpd};
    representation rep;

    solve_method sm, min_sm;

    int nrcols;
    int nrrows;
    int nrelem; // nrcols*nrrows;

    // the matrix elements:
    DL_Scalar *a;  // a[r][c] <=> a[r*nrcols+c]: doing our own address
               // calculation allows for optimisations.
    int asize; // the actual size of a (may be a bit larger than nrelem

    // lud representation:
    DL_Scalar *lu; // lower and upper matrices
    int *indx;     // permutation vector used by pivoting LUD-routines
    DL_Scalar d;   // number of row-exchanges in ludcmp odd or even;
    int bandw;     // bandwidth of the matrix (<=0: bandwidth not used)

    // limited (no sa-array) row indexed sparse storage representation:
    boolean *nonzero; // which element in the full representation is nonzero
    int *ijari;       // non-zero element start/stop indices in ijaci and ijami
    int *ijaci;       // non-zero element column indices
    int *ijami ;      // non-zero element matrix indices
    int nrnonzero;    // number of off-diagonal nonzero elements (length
                      // of ijaci and ijami

    // auxilary attributes for svd decomposition A=A' w V^T:
    DL_Scalar *u; // orthogonal nrrows x nrcols
    DL_Scalar *w; // diagonal nrcols x nrcols
    DL_Scalar *v; // square matrix nrcols x nrcols

    void  reptofull();
    void  full2riss();

    // without using the bandwidth
    int   ludcmp();
    void  lubksb(DL_largevector*, DL_largevector*);
    // using the bandwidth
    int   ludcmpbw();
    void  lubksbbw(DL_largevector*, DL_largevector*);

    boolean conjug_gradient(DL_largevector*, DL_largevector*);
    void  asolve(DL_largevector*,DL_largevector*);

    int   svdcmp();
    void  svbksb(DL_largevector*, DL_largevector*);

  public:

    inline int get_nrcols(){return nrcols;}
    inline int get_nrrows(){return nrrows;}
    
    void  assign(DL_largematrix*);
    void  assign(DL_matrix*);

    void  makezero();
    void  makeunit();

    void  resize(int,int);
    DL_Scalar  get(int,int);
    void  set(int,int,DL_Scalar);

    void  getsubmatrix(int,int,DL_largematrix*);
    void  setsubmatrix(int,int,DL_largematrix*);
    void  setsubmatrixnonzero(int, int, DL_largematrix*);
    void  setsubmatrixzero(int,int,int,int);
    void  setcolumn(int,DL_largevector*);
    void  setcolumn(int,DL_vector*);
    void  getcolumn(int,DL_largevector*);
    void  setrow(int,DL_largevector*);
    void  setrow(int,DL_vector*);
    void  getrow(int,DL_largevector*);

    void  plus(DL_largematrix*,DL_largematrix*);
    void  plusis(DL_largematrix*);
    void  minus(DL_largematrix*,DL_largematrix*);
    void  minusis(DL_largematrix*);

    void  times(DL_largematrix*,DL_largematrix*);
    void  times(DL_Scalar,DL_largematrix*);
    void  timesis(DL_Scalar);
    void  neg();
    void  times(DL_largevector*,DL_largevector*);
    void  transposetimes(DL_largevector*,DL_largevector*);

    int   get_bandwidth();
    int   get_nrnonzero();
  
    solve_method get_solve_method(){ return sm; };
    void  set_solve_method(solve_method);
    solve_method get_min_solve_method(){ return min_sm; };
    void  set_min_solve_method(solve_method);
  
    void  analyse_structure();
    boolean  prep_for_solve();
    boolean  solve(DL_largevector*, DL_largevector*);
  
    DL_Scalar det();

               DL_largematrix(int=0,int=0,solve_method=lud_bcksub); // constructor
               DL_largematrix(DL_largematrix*); // copy constructor
               ~DL_largematrix();               // destructor

    // for debugging: show myself:
    void show(void);
    void show_all(void);
};

inline DL_largematrix::DL_largematrix(int r, int c, solve_method _sm){
  nrrows=r;
  nrcols=c;
  asize=nrelem=r*c;
  if (asize<36) asize=36;
  a=new DL_Scalar[asize];
  rep=full;
  min_sm=lud_bcksub;
  sm=_sm;
  bandw=-1;
  d=1.0;
  indx=NULL;
  nonzero=NULL; ijari=ijaci=ijami=NULL;
  lu=u=w=v=NULL;
};

inline DL_largematrix::DL_largematrix(DL_largematrix *lm) {
  assign(lm);
}

inline DL_largematrix::~DL_largematrix() {
  if (a) delete[] a;
  if (lu) delete[] lu;
  if (indx) delete[] indx;
  if (nonzero) delete[] nonzero;
  if (ijari) {
    delete[] ijari;
    delete[] ijaci;
    delete[] ijami;
  }
  if (u) delete[] u;
  if (v) delete[] v;
  if (w) delete[] w;
}

inline void DL_largematrix::reptofull() {
  switch (rep) {
  case full: return;
  case riss:
    break;
  case lud:
  case ludb:
    bandw=-1;
    break;
  case svdcmpd:
    break;
  }        
  rep=full;
}

inline void DL_largematrix::makezero() {
  switch (rep) {
  case full:
    { for (int i=0;i<nrelem;i++) a[i]=0.0; }
    return;
  case riss:
    { int i;
      for (i=0; i<nrrows; i++) a[i*(nrcols+1)]=0.0;
      for (i=0;i<nrnonzero;i++) a[ijami[i]]=0.0;
    }
    return;
  case lud:
    rep=full;
    makezero();
    return;
  case ludb:
    rep=full;
    // only need to clear the band:
    {
      int r,c,ub,r_nrcols=0;
      for (r=0;r<nrrows;r++){
        ub=min(nrcols,r+2*bandw+1);
		for (c=max(0,r-2*bandw);c<ub;c++) a[r_nrcols+c]=0.0;
		r_nrcols+=nrcols;
      }
    }
    return;
  case svdcmpd:
    rep=full;
    makezero();
    return;
  }
}

inline void DL_largematrix::makeunit() {
  // pre: nrrows==nrcols;
  switch (rep) {
  case full:
  case riss:
    makezero();
    { int i; for (i=0; i<nrelem; i+=nrcols+1) a[i]=1.0; }
    return;
  case lud:
  case ludb:
    rep=full;
    makeunit();
    return;
  case svdcmpd:
    rep=full;
    makeunit();
  }
}

inline void DL_largematrix::resize(int r, int c) {
// This loses all info stored in the matrix!!!
  if (nonzero) delete[] nonzero; nonzero=NULL;
  bandw=-1;
  rep=full;
  nrrows=r;
  nrcols=c;
  nrelem=r*c;
  if ((nrelem>asize) || ((nrelem>36) && (2*nrelem<asize))){
    // enlarge `a' if the nr of elements won't fit
    // make `a' smaller if the size of `a' is larger than twice the
    // new number of elements, except for the cases where `a' is
    // already rather small (<=36 elements)
    delete[] a;
    asize=nrelem;
    if (asize==0) a=NULL;
    else a=new DL_Scalar[asize];
  }
  if (ijari){ delete[] ijari; ijari=NULL; }
  if (ijaci){ delete[] ijaci; ijaci=NULL; }
  if (ijami){ delete[] ijami; ijami=NULL; }
  if (lu) { delete[] lu; lu=NULL; }
  if (indx) { delete[] indx; indx=NULL; }
  if (u) { delete[] u; u=NULL; }
  if (v) { delete[] v; v=NULL; }
  if (w) { delete[] w; w=NULL; }
}

inline DL_Scalar DL_largematrix::get(int r,int c) {
// PRE: ((0<=r) && (r<nrrows) && (0<=c) && (c<nrcols))
//      && (rep==full/riss)
    return a[r*nrcols+c];
}

inline void DL_largematrix::set(int r,int c,DL_Scalar f) {
// PRE: ((0<=r) && (r<nrrows) && (0<=c) && (c<nrcols))
//      && (rep==full/riss)
    a[r*nrcols+c]=f;
}

inline void DL_largematrix::getsubmatrix(int r,int c,DL_largematrix* lm) {
// PRE: lm && (r+lm->nrrows<=nrrows) && (c+lm->nrcols<=nrcols)
//      && (rep==full/riss)
  lm->reptofull();
  int ri, ci, ri_lm_nrcols=0;   // INV: ri_lm_nrcols==ri*lm->nrcols
  int r_ri_nrcols_c=r*nrcols+c; // INV: r_ri_nrcols_c==(r+ri)*nrcols+c
  for (ri=0; ri<lm->nrrows; ri++) {
    for (ci=0; ci<lm->nrcols; ci++)
      lm->a[ri_lm_nrcols+ci]=a[r_ri_nrcols_c+ci];
    ri_lm_nrcols+=lm->nrcols;
    r_ri_nrcols_c+=nrcols;
  }
}

inline void DL_largematrix::setsubmatrix(int r,int c,DL_largematrix* lm){
// PRE: lm && (r+lm->nrrows<=nrrows) && (c+lm->nrcols<=nrcols)
//      && (rep==full/riss) && (lm->rep==full/riss)
  int ri, ci, ri_lm_nrcols=0;   // INV: ri_lm_nrcols==ri*lm->nrcols
  int r_ri_nrcols_c=r*nrcols+c; // INV: r_ri_nrcols_c==(r+ri)*nrcols+c
  for (ri=0; ri<lm->nrrows; ri++) {
    for (ci=0; ci<lm->nrcols; ci++)
      a[r_ri_nrcols_c+ci]=lm->a[ri_lm_nrcols+ci];
    ri_lm_nrcols+=lm->nrcols;
    r_ri_nrcols_c+=nrcols;
  }
    
}

inline void DL_largematrix::setsubmatrixzero(int r,int c,int i,int j){
// PRE: lm && (r+i<=nrrows) && (c+j<=nrcols)
  switch (rep) {
  case full: break;
  case riss:
    delete[] ijari; ijari=NULL;
    delete[] ijaci; ijaci=NULL;
    delete[] ijami; ijami=NULL;
    rep=full;
    break;
  case lud:
  case ludb:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::setsubmatrixzero:\n Can not set elements of a LU decomposed matrix\n");
    return;
  case svdcmpd:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::setsubmatrixzero:\n Can not set elements of an SV decomposed matrix\n");
    return;
  }
  int ri,ci,ri_nrcols; // INV: ri_nrcols==ri*nrcols
  if (!nonzero)
    if (nrrows==nrcols) {
      nonzero=new boolean[nrelem];
      for (ri=0;ri<nrelem;ri++) nonzero[ri]=FALSE;
    }
  if (nonzero)
    for(ri=r; ri<r+i; ri++) {
      ri_nrcols=ri*nrcols;
      for (ci=c; ci<c+j; ci++)
	nonzero[ri_nrcols+ci]=FALSE;
    }
  for(ri=r; ri<r+i; ri++) {
    ri_nrcols=ri*nrcols;
    for (ci=c; ci<c+j; ci++)
      a[ri_nrcols+ci]=0.0;
  }
}

inline void DL_largematrix::setsubmatrixnonzero(int r,int c,DL_largematrix* lm){
// PRE: lm && (r+lm->nrrows<=nrrows) && (c+lm->nrcols<=nrcols)
//      && (lm->rep!=lud/b)
  switch (rep) {
  case full: break;
  case riss:
    delete[] ijari; ijari=NULL;
    delete[] ijaci; ijaci=NULL;
    delete[] ijami; ijami=NULL;
    rep=full;
    break;
  case lud:
  case ludb:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::setsubmatrixnonzero:\n Can not set elements of a LU decomposed matrix\n");
    return;
  case svdcmpd:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::setsubmatrixnonzero:\n Can not set elements of an SV decomposed matrix\n");
    return;
  }
  int ri,ci;
  if (!nonzero)
    if (nrrows==nrcols) {
      nonzero=new boolean[nrelem];
      for (ri=0;ri<nrelem;ri++) nonzero[ri]=FALSE;
    }
    else
      DL_dsystem->get_companion()->Msg("DL_largematrix::setsubmatrixnonzero is only effective for square matrices. Use DL_largematrix::setsubmatrix instead\n");
  
  register int idx0;
  register int idx1;
  if (nonzero) {
    idx0=r*nrcols; // INV: idx0==ri*nrcols;
    for(ri=r; ri<r+lm->nrrows; ri++) {
      for (ci=c; ci<c+lm->nrcols; ci++)
	nonzero[idx0+ci]=TRUE;
      idx0+=nrcols;
    }
  }
  idx0=r*nrcols+c; // INV: idx0==(r+ri)*nrcols+c
  idx1=0;          // INV: idx1==ri*lm->nrcols;
  for (ri=0; ri<lm->nrrows; ri++) {
    for (ci=0; ci<lm->nrcols; ci++)
      a[idx0+ci]=lm->a[idx1+ci];
    idx0+=nrcols;
    idx1+=lm->nrcols;
  }
}

inline void DL_largematrix::setcolumn(int c,DL_largevector* lv) {
// PRE: lv->dim>=nrrows && 0<=c<nrcols
  switch (rep) {
  case full:
  case riss:
    {
      register int r;
      register int r_nrcols_c=c;  // INV: r_nrcols_c==r*nrcols+c
      for(r=0;r<nrrows;r++) {
	a[r_nrcols_c]=lv->get(r);
	r_nrcols_c+=nrcols;
      }
    }
    return;
  case lud:
  case ludb:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::setcolumn:\n Can not set elements of a LU decomposed matrix\n");
    return;
  case svdcmpd:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::setcolumn:\n Can not set elements of an SV decomposed matrix\n");
    return;
  }
}

inline void DL_largematrix::setcolumn(int c,DL_vector* vv) {
// PRE: 3<=nrrows && 0<=c<nrcols
//      && (rep==full/riss)
  register int r_nrcols_c=c;  // INV: r_nrcols_c==r*nrcols+c
  a[r_nrcols_c]=vv->x; r_nrcols_c+=nrcols;
  a[r_nrcols_c]=vv->y; r_nrcols_c+=nrcols;
  a[r_nrcols_c]=vv->z;
}

inline void DL_largematrix::getcolumn(int c,DL_largevector* lv) {
// PRE: lv->dim>=nrrows && 0<=c<nrcols
  register int r;
  register int r_nrcols_c=c;  // INV: r_nrcols_c==r*nrcols+c
  for(r=0;r<nrrows;r++) {
    lv->set(r,a[r_nrcols_c]);
    r_nrcols_c+=nrcols;
  }
}

inline void DL_largematrix::setrow(int r,DL_largevector* lv){
// PRE: lv->dim>=nrcols && 0<=r<nrrows
  switch (rep) {
  case full:
    {
      register int c,r_nrcols=r*nrcols;
      for(c=0;c<nrcols;c++) a[r_nrcols+c]=lv->get(c);
    }
    return;
  case riss:
    {
      a[r*(nrcols+1)]=lv->get(r);
      for(int i=ijari[r]; i<ijari[r+1]; i++)
	a[ijami[i]]=lv->get(ijaci[i]);
    }
    return;
  case lud:
  case ludb:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::setrow:\n Can not set elements of a LU decomposed matrix\n");
    return;
  case svdcmpd:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::setrow:\n Can not set elements of an SV decomposed matrix\n");
    return;
  }
}

inline void DL_largematrix::setrow(int r,DL_vector* vv){
// PRE: 3<=nrcols && 0<=r<nrrows
//      && (rep==full/riss)
  register int r_nrcols=r*nrcols;
  a[r_nrcols++]=vv->x;
  a[r_nrcols++]=vv->y;
  a[r_nrcols]=vv->z;
}

inline void DL_largematrix::getrow(int r,DL_largevector* lv){
// PRE: lv->dim>=nrcols && 0<=r<nrrows
  register int c,r_nrcols=r*nrcols;
  for(c=0;c<nrcols;c++) lv->set(c,a[r_nrcols+c]);
}


inline void DL_largematrix::plus(DL_largematrix *lm, DL_largematrix *nlm) {
// PRE: lm && nlm &&
//      (nrrows==lm->nrrows==nlm->nrrows) && (nrcols==lm->nrcols==nlm->nrcols)
//      && (rep==full/riss) && (lm->rep==full/riss)
  nlm->reptofull();
  for (int i=0; i<nrelem; i++) nlm->a[i]=a[i]+lm->a[i];
}

inline void DL_largematrix::plusis(DL_largematrix *lm) {
// PRE: lm &&
//      (nrrows==lm->nrrows) && (nrcols==lm->nrcols)
//      && (rep==full/riss) && (lm->rep==full/riss)
  reptofull();
  for (int i=0; i<nrelem; i++) a[i]+=lm->a[i];
}

inline void DL_largematrix::minus(DL_largematrix *lm, DL_largematrix *nlm) {
// PRE: lm && nlm &&
//      (nrrows==lm->nrrows==nlm->nrrows) && (nrcols==lm->nrcols==nlm->nrcols)
//      && (rep==full/riss) && (lm->rep==full/riss)
  nlm->reptofull();
  for (int i=0; i<nrelem; i++) nlm->a[i]=a[i]-lm->a[i];
}

inline void DL_largematrix::minusis(DL_largematrix *lm) {
// PRE: lm &&
//      (nrrows==lm->nrrows) && (nrcols==lm->nrcols)
//      && (rep==full/riss) && (lm->rep==full/riss)
  reptofull();
  for (int i=0; i<nrelem; i++) a[i]-=lm->a[i];
}

inline void DL_largematrix::times(DL_largematrix *lm, DL_largematrix *nlm) {
//PRE: lm && nlm && (nrrows==nlm->nrrows) && (lm->nrcols==nlm->nrcols)
//     && (nrcols==lm->nrrows)
//     && (rep==full/riss) && (lm->rep==full/riss)
  register int ri,ci,i;
  register int ri_nrcols=0;
  int i_lm_nrcols_ci=0;
  int ri_nlm_nrcols=0;
  DL_Scalar temp;
  nlm->reptofull();
  for (ri=0; ri<nlm->nrrows; ri++) {
    for (ci=0; ci<nlm->nrcols; ci++) {
      temp=0.0; i_lm_nrcols_ci=ci;
      for (i=0; i<nrcols; i++) {
	temp+=a[ri_nrcols+i]*lm->a[i_lm_nrcols_ci];
	i_lm_nrcols_ci+=lm->nrcols;
      }
      nlm->a[ri_nlm_nrcols+ci]=temp;
    }
    ri_nrcols+=nrcols;
    ri_nlm_nrcols+=nlm->nrcols;
  }
}

inline void DL_largematrix::times(DL_Scalar f, DL_largematrix *nlm) {
// PRE: nlm && (nrrows==nlm->nrrows) && (nrcols==nlm->nrcols)
  nlm->reptofull();
  for(int i=0; i<nrelem; i++) nlm->a[i]=f*a[i];
}

inline void DL_largematrix::timesis(DL_Scalar f) {
  switch (rep) {
  case full:
    { for(int i=0;i<nrelem; i++) a[i]*=f; }
    return;
  case riss:
    {
      int r,i;
      for (r=0;r<nrrows;r++) {
	a[r*(nrcols+1)]*=f;
	for (i=ijari[r];i<ijari[r+1];i++)
	  a[ijami[i]]*=f;
      }
    }
    return;
  case lud:
  case ludb:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::timesis:\n Can not multiply a LU decomposed matrix\n");
    return;
  case svdcmpd:
    { for (int i=0;i<nrcols;i++) w[i]*=f; }
    return;
  }
}

inline void DL_largematrix::neg() {
  switch (rep) {
  case full:
    { for(int i=0; i<nrelem; i++) a[i]=-a[i]; }
    return;
  case riss:
    {
      int r,i,k;
      for (r=0;r<nrrows;r++) {
	k=r*(nrcols+1);
	a[k]=-a[k];
	for (i=ijari[r];i<ijari[r+1];i++) {
	  k=ijami[i];
	  a[k]=-a[k];
	}
      }
    }
    return;
  case lud:
  case ludb:
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::neg:\n Can not negate a LU decomposed matrix\n");
    return;
  case svdcmpd:
    { for (int i=0;i<nrcols;i++) w[i]=-w[i]; }
    return;
  }
}

inline void DL_largematrix::times(DL_largevector *lv,DL_largevector *nlv) {
//PRE: lv && nlv && (nrrows==lv->dim) && (nrcols==nlv->dim)
  if (rep==riss) {
    int r,i;
      DL_Scalar temp;
      for (r=0;r<nrrows;r++) {
	temp=a[r*(nrcols+1)]*lv->get(r);
	for (i=ijari[r]; i<ijari[r+1]; i++)
	  temp+=a[ijami[i]]*lv->get(ijaci[i]);
	nlv->set(r,temp);
      }
    }
  else {
    int r,c,r_nrcols;
    DL_Scalar temp;
    for (r=0; r<nrrows; r++) {
      temp=0.0; r_nrcols=r*nrcols;
      for (c=0; c<nrcols; c++) temp+=a[r_nrcols+c]*lv->get(c);
      nlv->set(r,temp);
    }
  }
}

inline void DL_largematrix::transposetimes(DL_largevector *lv,DL_largevector *nlv) {
//PRE: lv && nlv && (nrrows==lv->dim) && (nrcols==nlv->dim)
  if (rep==riss) {
    int r,i;
    for (r=0;r<nrrows;r++) nlv->set(r,a[r*(nrcols+1)]*lv->get(r));
    for (r=0;r<nrrows;r++) {
      for (i=ijari[r]; i<ijari[r+1]; i++) {
	nlv->v[ijaci[i]]+=a[ijami[i]]*lv->get(r);
	// not so nice to use nlv's implementation, but += is faster ...
      }
    }
  }
  else {
    int r,c,c_nrcols_r;
    DL_Scalar temp;
    for (r=0; r<nrrows; r++) {
      temp=0.0; c_nrcols_r=r;
      for (c=0; c<nrcols; c++) {
	temp+=a[c_nrcols_r]*lv->get(c);
	c_nrcols_r+=nrcols;
      }
      nlv->set(r,temp);
    }
  }
}

inline void DL_largematrix::asolve(DL_largevector *b, DL_largevector *x){
// auxilary method for conjug_grad
// PRE: nrrows==nrcols==b->dim==x->dim && rep==full/riss
  for (int i=0; i<nrrows ; i++) {
    DL_Scalar tmp=a[i*(nrcols+1)];
    x->set(i,(tmp!=0.0 ? b->get(i)/tmp: b->get(i)));
  }
}

inline void DL_largematrix::assign(DL_matrix* m) {
// PRE: m && (nrcols==nrrows==3) && rep==full/riss
  a[0]=m->c0.x; a[1]=m->c1.x; a[2]=m->c2.x;
  a[3]=m->c0.y; a[4]=m->c1.y; a[5]=m->c2.y;
  a[6]=m->c0.z; a[7]=m->c1.z; a[8]=m->c2.z;
}

inline int DL_largematrix::get_bandwidth(){
  if (bandw>=0) return bandw;
  if (!nonzero) {
    if (nrelem>0)
      DL_dsystem->get_companion()->Msg("Error: zero-structure not known in get_bandwidth()\n");
    return nrcols;
  }
  bandw=0;
  int r,c,r_nrcols=0;
  for (r=0;r<nrrows;r++) {
    for (c=0;c<nrcols;c++)
      if (nonzero[r_nrcols+c])
	bandw=max(bandw,abs(r-c));
    r_nrcols+=nrcols;
  }
  return bandw;
}

inline int DL_largematrix::get_nrnonzero(){
  if (!nonzero) {
    if (nrelem>0)
      DL_dsystem->get_companion()->Msg("Error: zero-structure not known in calc_nrnonzero()\n");
    nrnonzero=nrelem;
    return nrnonzero;
  }

  int i,j,i_nrcols=0;
  
  nrnonzero=0;
  for (i=0; i<nrrows; i++) {
    for (j=0; j<nrrows; j++)
      if ((i!=j) && nonzero[i_nrcols+j]) nrnonzero++;
    i_nrcols+=nrcols;
  }
  return nrnonzero;
}

#endif
