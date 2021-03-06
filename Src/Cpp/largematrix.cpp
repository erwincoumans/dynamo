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
// filename	: largematrix.cpp
// description	: non-inline methods of class DL_largematrix
// author	: Bart Barenbrug     Sept 1997
//

#include "largematrix.h"
//#define DEBUG

// use pivoting in LU decomposition/backward substitution or not:
//#define PIVOT

// ************************** //
// non-inline member fuctions //
// ************************** //

void DL_largematrix::assign(DL_largematrix *lm) {
// PRE: lm
  resize(lm->nrrows,lm->nrcols);
  int i;
  for (i=0; i<nrelem; i++) a[i]=lm->a[i];
  if (lm->nonzero) {
    if (nonzero) delete[] nonzero;
    nonzero=new boolean[nrelem];
    for (i=0; i<nrelem; i++) nonzero[i]=lm->nonzero[i];
  }
  nrnonzero=lm->nrnonzero;
  bandw=lm->bandw;
  d=lm->d;
  rep=lm->rep;
  sm=lm->sm;
  min_sm=lm->min_sm;
  switch (lm->rep) {
  case full: break;
  case riss:
    rep=riss;
    if (!ijari) ijari=new int[nrrows+1];
    if (!ijaci) ijaci=new int[nrnonzero];
    if (!ijami) ijami=new int[nrnonzero];
    for (i=0; i<=nrrows; i++) ijari[i]=lm->ijari[i];
    for (i=0; i<nrnonzero; i++) ijaci[i]=lm->ijaci[i];
    for (i=0; i<nrnonzero; i++) ijami[i]=lm->ijami[i];
    break;
  case lud:
    rep=lud;
    if (!lu) lu=new DL_Scalar[nrelem];
    for (i=0;i<nrelem;i++) lu[i]=lm->lu[i];
    if (!indx) indx=new int[nrrows];
    for (i=0;i<nrrows;i++) indx[i]=lm->indx[i];
    break;
  case ludb:
    rep=ludb;
    if (!lu) lu=new DL_Scalar[nrelem];
    for (i=0;i<nrelem;i++) lu[i]=lm->lu[i];
    if (!indx) indx=new int[nrrows];
    for (i=0;i<nrrows;i++) indx[i]=lm->indx[i];
    break;
  case svdcmpd:
    rep=svdcmpd;
    if (!u) u=new DL_Scalar[nrelem];
    if (!w) w=new DL_Scalar[nrcols];
    if (!v) v=new DL_Scalar[nrcols*nrcols];
    for (i=0;i<nrelem;i++) u[i]=lm->u[i];
    for (i=0;i<nrcols;i++) w[i]=lm->w[i];
    for (i=0;i<nrcols*nrcols;i++) v[i]=lm->v[i];
    break;
  }
}

void DL_largematrix::full2riss() {
// PRE: rep==full
  if (!nonzero) {
    DL_dsystem->get_companion()->Msg("Error: zero-structure not known in full2riss()\n");
    return;
  }

  int i,j,k,i_nrcols=0;
  
  if (!ijari) ijari=new int[nrrows+1];
  if (!ijaci) ijaci=new int[nrnonzero];
  if (!ijami) ijami=new int[nrnonzero];

  ijari[0]=k=i_nrcols=0;
  for(i=0;i<nrrows;i++) {
    for(j=0;j<nrrows;j++) {
      if (nonzero[i_nrcols+j] && (i!=j)) {
	ijaci[k]=j;
	ijami[k]=i_nrcols+j;
	k++;
      }
    }
    ijari[i+1]=k;
    i_nrcols+=nrcols;
  }

  rep=riss;
}

#ifdef PIVOT

#define TINY 1.0e-10
int DL_largematrix::ludcmp(){
// Calculates the LU-decomposion of the matrix. The lower and
// upper matrices are stored in a, while the permutation vector
// is stored in indx, and d stores if there have been an odd or
// even number of row-exchanges (used in determinant determination)

// If the decomposition went ok, a negative integer is returned.
// Otherwise the index of the offending row is returned. In the latter
// case the contents of the matrix will be undefined

// PRE: nrrows==nrcols && rep==full

  static DL_largevector vv;
  int i, imax, j, k, nrcols_i, nrcols_j_j=0;
  DL_Scalar big,dum,sum,temp;

  if (!lu) lu=new DL_Scalar[nrelem];
  if (!indx) indx=new int[nrcols];
  vv.resize(nrcols);
  d=1.0;

  nrcols_i=0;
  for(i=0;i<nrcols;i++) {
    big=0.0;
    for(j=0;j<nrcols;j++)
      if ((temp=fabs(a[nrcols_i+j])) > big) big=temp;
    if (big<TINY) {
      delete[] indx; indx=NULL;
      return i;
    }
    vv.set(i,1.0/big);
    nrcols_i+=nrcols;
  }
  for(j=0;j<nrcols;j++) {
      nrcols_i=0;
      for(i=0;i<j;i++) {
          sum=a[nrcols_i+j];
          int nrcols_k_j=j;   // INV: nrcols_k_j==nrcols*k+j
          for(k=0;k<i;k++) {
              sum-=lu[nrcols_i+k]*lu[nrcols_k_j];
	      nrcols_k_j+=nrcols;
          }
          lu[nrcols_i+j]=sum;
          nrcols_i+=nrcols;
      }
      big=0.0;
      for(i=j;i<nrcols;i++) {
          int nrcols_k_j=j;   // INV: nrcols_k_j==nrcols*k+j
          sum=a[nrcols_i+j];
          for(k=0;k<j;k++) {
              sum-=lu[nrcols_i+k]*lu[nrcols_k_j];
	      nrcols_k_j+=nrcols;
          }
          lu[nrcols_i+j]=sum;
          if ( (dum=vv.get(i)*fabs(sum)) >=big) {
	      big=dum;
	      imax=i;
	  }
	  nrcols_i+=nrcols;
      }
      if (j != imax) {
          int nrcols_imax_k=nrcols*imax; // INV: nrcols_imax_k==nrcols*imax+k
          int nrcols_j_k=nrcols*j;       // INV: nrcols_j_k==nrcols*j+k
          for(k=0;k<nrcols;k++) {
	      dum=lu[nrcols_imax_k];
	      lu[nrcols_imax_k++]=lu[nrcols_j_k];
	      lu[nrcols_j_k++]=dum;
	  }
	  d=-d;
	  vv.set(imax,vv.get(j));
      }
      indx[j]=imax;
      if (lu[nrcols_j_j] == 0.0) lu[nrcols_j_j]=TINY;
      if (j!=nrcols-1) {
          dum=1.0/(lu[nrcols_j_j]);
	  int nrcols_i_j=nrcols*(j+1)+j;
	  for (i=j+1;i<nrcols;i++) {
	    lu[nrcols_i_j]*=dum;
	    nrcols_i_j+=nrcols;
	  }
      }
      nrcols_j_j+=nrcols+1;
  }
  rep=lud;
  return -1;
}


int DL_largematrix::ludcmpbw(){
// Calculates the LU-decomposion of the matrix using the bandwidth.
// The lower and upper matrices are stored in a, while the permutation vector
// is stored in indx, and d stores if there have been an odd or
// even number of row-exchanges (used in determinant determination)

// PRE: nrrows==nrcols && rep==full

  static DL_largevector vv;
  DL_Scalar sum,big,dum,temp;
  int i,j,imax,k,lb,ub,nrcols_i=0; // INV: nrcols_i==nrcols*i

  if (!lu) lu=new DL_Scalar[nrelem];
  for(i=0;i<nrelem;i++) lu[i]=0.0;
  if (!indx) indx=new int[nrcols];
  vv.resize(nrcols);
  d=1.0;
  
  for(i=0;i<nrcols;i++) {
    big=0.0;
    lb=max(0,i-bandw);
    ub=min(nrcols,i+bandw+1);
    for(j=lb;j<ub;j++) {
      if ((temp=fabs(a[nrcols_i+j])) > big) big=temp;
    }
    if (big<TINY) {
      delete[] indx; indx=NULL;
      return i;
    }
    vv.set(i,1.0/big);
    nrcols_i+=nrcols;
  }
  int nrcols_j_j=0; // INV: nrcols_j_j==nrcols*j+j
  for(j=0;j<nrcols;j++) {
    lb=max(0,j-2*bandw-1);
    nrcols_i=lb*nrcols;
    for(i=lb;i<j;i++) {
      int nrcols_k_j=j;   // INV: nrcols_k_j==nrcols*k+j
      sum=a[nrcols_i+j];
      for(k=0;k<i;k++) {
	sum-=lu[nrcols_i+k]*lu[nrcols_k_j];
	nrcols_k_j+=nrcols;
      }
      lu[nrcols_i+j]=sum;
      nrcols_i+=nrcols;
    }
    big=0.0;
    ub=min(nrcols,j+bandw+1);
    nrcols_i=j*nrcols;
    for(i=j;i<ub;i++) {
      int nrcols_k_j=j;   // INV: nrcols_k_j==nrcols*k+j
      sum=a[nrcols_i+j];
      for(k=0;k<j;k++) {
	sum-=lu[nrcols_i+k]*lu[nrcols_k_j];
	nrcols_k_j+=nrcols;
      }
      lu[nrcols_i+j]=sum;
      if ( (dum=vv.get(i)*fabs(sum)) >=big) {
	big=dum;
	imax=i;
      }
      nrcols_i+=nrcols;
    }
    if (j != imax) {
      int nrcols_imax_k=nrcols*imax; // INV: nrcols_imax_k==nrcols*imax+k
      int nrcols_j_k=nrcols*j;       // INV: nrcols_j_k==nrcols*j+k
      for(k=0;k<nrcols;k++) {
	dum=a[nrcols_imax_k];
	lu[nrcols_imax_k]=lu[nrcols_j_k];
	lu[nrcols_j_k]=dum;
	nrcols_imax_k++;
	nrcols_j_k++;
      }
      d=-d;
      vv.set(imax,vv.get(j));
    }
    indx[j]=imax;
    if (lu[nrcols_j_j] == 0.0) lu[nrcols_j_j]=TINY;
    if (j+1!=nrcols) {
      dum=1.0/(lu[nrcols_j_j]);
      int nrcols_i_j=(j+1)*nrcols+j;
      for (i=j+1;i<nrcols;i++) {
	lu[nrcols_i_j]*=dum;
	nrcols_i_j+=nrcols;
      }
    }
    nrcols_j_j+=nrcols+1;
  }
  rep=ludb;
  return -1;
}
#undef TINY

void DL_largematrix::lubksb(DL_largevector *x, DL_largevector *b){
// Using the LU decomposition stored in this matrix, solves x from
// self x=b.

// PRE: nrrows==nrcols==b->dim && rep==ludb

  int i, ii=nrcols, ip, j;
  DL_Scalar sum;

  int nrcols_i=0; // INV: nrcols_i==nrcols*i;
  for(i=0;i<nrcols;i++) {
    ip=indx[i];
    sum=b->get(ip);
    x->set(ip,b->get(i));
    if (ii!=nrcols)
      for(j=ii;j<i;j++)
	sum-=lu[nrcols_i+j]*x->get(j);
    else if (sum) ii=i;
    x->set(i,sum);
    nrcols_i+=nrcols;
  }
  
  for (i=nrcols-1;i>=0;i--) {
    nrcols_i-=nrcols;
    sum=x->get(i);
    for(j=i+1;j<nrcols;j++)
      sum-=lu[nrcols_i+j]*x->get(j);
    x->set(i,sum/a[nrcols_i+i]);
  }
}

void DL_largematrix::lubksbbw(DL_largevector *x,DL_largevector *b){
// Using the LU decomposition stored in this matrix, solves x from
// self x=b.

// PRE: nrrows==nrcols==b->dim  && rep==ludb

  DL_Scalar sum;
  int i,j,ip,ub,ii=nrcols,nrcols_i=0; // INV: nrcols_i==nrcols*i;

  for(i=0;i<nrcols;i++) {
    ip=indx[i];
    sum=b->get(ip);
    x->set(ip,b->get(i));
    if (ii!=nrcols)
      for(j=ii;j<i;j++)
	sum-=lu[nrcols_i+j]*x->get(j);
    else if (sum) ii=i;
    x->set(i,sum);
    nrcols_i+=nrcols;
  }
  for (i=nrcols-1;i>=0;i--) {
    nrcols_i-=nrcols;
    sum=x->get(i);
//    ub=min(nrcols,i+2*bandw+1);
    ub=min(nrcols,i+bandw+1);
    for(j=i+1;j<ub;j++) sum-=lu[nrcols_i+j]*x->get(j);
    x->set(i,sum/lu[nrcols_i+i]);
  }
}

#else  // notdef PIVOT

#define TINY 1.0e-10
int DL_largematrix::ludcmp(){
// Calculates the LU-decomposion of the matrix. The lower and
// upper matrices are stored in lu

// If the decomposition went ok, a negative integer is returned.
// Otherwise the index of the offending row is returned. In the latter
// case the contents of the matrix will be undefined

// PRE: nrrows==nrcols && rep==full

  int i, j, k;
  DL_Scalar sum,dum;
  if (!lu) lu=new DL_Scalar[nrelem];

  int nrcols_i=0;   // INV: nrcols_i==nrcols*i
  int nrcols_j_j=0; // INV: nrcols_j_j==nrcols*j+j
  
  for(j=0;j<nrcols;j++) {
    nrcols_i=0;
    for(i=0;i<=j;i++) {
      int nrcols_k_j=j;   // INV: nrcols_k_j==nrcols*k+j
      sum=a[nrcols_i+j];
      for(k=0;k<i;k++) {
	sum-=lu[nrcols_i+k]*lu[nrcols_k_j];
	nrcols_k_j+=nrcols;
      }
      lu[nrcols_i+j]=sum;
      nrcols_i+=nrcols;
    }      
    if (fabs(lu[nrcols_j_j]) < TINY) return j;
    dum=1.0/lu[nrcols_j_j];
    for (i=j+1;i<nrcols;i++) {
      int nrcols_k_j=j;
      sum=a[nrcols_i+j];
      for (k=0;k<j;k++) {
	sum-=lu[nrcols_i+k]*lu[nrcols_k_j];
	nrcols_k_j+=nrcols;
      }
      lu[nrcols_i+j]=sum*dum;
      nrcols_i+=nrcols;
    }
    nrcols_j_j+=nrcols+1;
  }
  rep=lud;
  return -1;
}

int DL_largematrix::ludcmpbw(){
// Calculates the LU-decomposion of the matrix using the bandwidth.
// The lower and upper matrices are stored in a.

// PRE: nrrows==nrcols && rep==full

  int i, j, k, lb, ub;
  DL_Scalar sum,dum;
  if (!lu) lu=new DL_Scalar[nrelem];
  for(i=0;i<nrelem;i++) lu[i]=0.0;

  int nrcols_i=0;   // INV: nrcols_i==nrcols*i
  int nrcols_j_j=0; // INV: nrcols_j_j==nrcols*j+j
  
  for(j=0;j<nrcols;j++) {
    lb=max(0,j-bandw);
    nrcols_i=lb*nrcols;
    for(i=lb;i<=j;i++) {
      int nrcols_k_j=lb*nrcols+j;   // INV: nrcols_k_j==nrcols*k+j
      sum=a[nrcols_i+j];
      for(k=lb;k<i;k++) {
	sum-=lu[nrcols_i+k]*lu[nrcols_k_j];
	nrcols_k_j+=nrcols;
      }
      lu[nrcols_i+j]=sum;
      nrcols_i+=nrcols;
    }      
    if (fabs(lu[nrcols_j_j]) < TINY) return j;
    dum=1.0/lu[nrcols_j_j];
    ub=min(nrcols,j+bandw+1);
    for (i=j+1;i<ub;i++) {
      int nrcols_k_j=lb*nrcols+j;
      sum=a[nrcols_i+j];
      for (k=lb;k<j;k++) {
	sum-=lu[nrcols_i+k]*lu[nrcols_k_j];
	nrcols_k_j+=nrcols;
      }
      lu[nrcols_i+j]=sum*dum;
      nrcols_i+=nrcols;
    }
    nrcols_j_j+=nrcols+1;
  }
  rep=ludb;
  return -1;
}

#undef TINY

void DL_largematrix::lubksb(DL_largevector *x, DL_largevector *b){
// Using the LU decomposition stored in this matrix, solves x from
// self x=b.

// PRE: nrrows==nrcols==b->dim && rep==ludb

  DL_Scalar sum;
  int i,j,nrcols_i=nrcols; // INV: nrcols_i == nrcols*i

  x->set(0,b->get(0));
  for(i=1;i<nrcols;i++) {
    sum=b->get(i);
    for(j=0;j<i;j++) sum-=lu[nrcols_i+j]*x->get(j);
    x->set(i,sum);
    nrcols_i+=nrcols;
  }
  
  for (i=nrcols-1;i>=0;i--) {
    nrcols_i-=nrcols;
    sum=x->get(i);
    for(j=i+1;j<nrcols;j++) sum-=lu[nrcols_i+j]*x->get(j);
    x->set(i,((DL_Scalar)sum)/lu[nrcols_i+i]);
  }
}

void DL_largematrix::lubksbbw(DL_largevector *x, DL_largevector *b){
// Using the LU decomposition stored in this matrix, solves x from
// self x=b. x is returned in b.

// PRE: nrrows==nrcols==b->dim  && rep==ludb

  int i, j, ub;
  DL_Scalar sum;

  int nrcols_i=0; // INV: nrcols_i == nrcols*i

  for(i=0;i<nrcols;i++) {
    sum=b->get(i);
    for(j=max(0,i-bandw);j<i;j++) sum-=lu[nrcols_i+j]*x->get(j);
    x->set(i,sum);
    nrcols_i+=nrcols;
  }
  
  for (i=nrcols-1;i>=0;i--) {
    nrcols_i-=nrcols;
    sum=x->get(i);
    ub=min(nrcols,i+bandw+1);
    for(j=i+1;j<ub;j++) sum-=lu[nrcols_i+j]*x->get(j);
    x->set(i,sum/lu[nrcols_i+i]);
  }
}
  
#endif //def PIVOT

DL_Scalar DL_largematrix::det() {
// Calculates the determinant of the matrix using the LU decomposition.
//PRE: nrcols=nrrows
  DL_Scalar det=d;
  if ((rep!=lud) && (rep!=ludb)) {
    representation org_rep=rep;
    ludcmp();
    rep=org_rep;
  }
  for (int i=0; i<nrcols; i++) det*=lu[i*(nrcols+1)];
  return det;
}

#define TOL 0.1
boolean DL_largematrix::conjug_gradient(DL_largevector *x, DL_largevector *b){
// Solves Ax=b using conjugate gradient
// returns if the solution was diverging |Ax-b|>|b|
  static DL_largevector p,pp, r,rr, z,zz;

  if ((rep==lud) || (rep==ludb)) {
    DL_dsystem->get_companion()->Msg("DL_largematrix::conjug_grad not implemented for LU decomposed matrices\n");
    return FALSE;
  }

  int j, iter=0, itmax=nrrows;
  DL_Scalar ak, akden, bk, bkden=0, bknum, bnrm;

  p.resize(nrrows);
  pp.resize(nrrows);
  r.assign(b);
  rr.assign(b);
  z.resize(nrrows);
  zz.resize(nrrows);

  bnrm=b->norm();
  x->makezero();
  asolve(&r,&z);
  while (iter<itmax) {
    ++iter;
    asolve(&rr,&zz);
    bknum=z.inprod(&rr);
    if (iter==1) {
      p.assign(&z);
      pp.assign(&zz);
    }
    else {
      bk=bknum/bkden;
      for (j=0;j<nrrows;j++) {
	p.set(j,bk*p.get(j)+z.get(j));
	pp.set(j,bk*pp.get(j)+zz.get(j));
      }
    }
    bkden=bknum;
    times(&p,&z);
    akden=z.inprod(&pp);
    ak=bknum/akden;
    transposetimes(&pp,&zz);
    for(j=0;j<nrrows;j++) { // not very nice to use the implementation
      // of x, r and rr here, but += and -= are faster...
      x->v[j]+=ak*p.get(j);
       r.v[j]-=ak*z.get(j);
      rr.v[j]-=ak*zz.get(j);
    }
    asolve(&r,&z);

    if (r.norm()<=(TOL*bnrm)) break;
  }
  return (r.norm()>bnrm);
}

static DL_Scalar at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
int DL_largematrix::svdcmp(){
// calculates an singular value decompisition of this matrix:
// A=U w V^T
// where U is nrrows x nrcols and has orthonormal columns (stored in a)
// where w is nrcols x nrcols diagonal with the singular values
// where V is nrcols x nrcols and has orthonormal columns
// pre: nrrows>=nrcols && rep==full
// returns te number of singular values after postprocessing w to
// nullify too small elements in w
  int flag,i,its,j,jj,k,l,nm;
  DL_Scalar c,f,h,s,x,y,z;
  DL_Scalar anorm=0.0,g=0.0,scale=0.0;
  DL_Scalar *rv1;

  if (nrrows < nrcols) {
    DL_dsystem->get_companion()->Msg("Error: DL_largematrix::svdcmp: more columns than rows!\n");
    return 1;
  }

  rep=svdcmpd;
  if (!u) u=new DL_Scalar[nrelem];
  if (!w) w=new DL_Scalar[nrcols];
  if (!v) v=new DL_Scalar[nrcols*nrcols];
  rv1=new DL_Scalar[nrcols];

  for (i=0;i<nrelem;i++) u[i]=a[i];
  for (i=0;i<nrcols;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < nrrows) {
      register int k_nrcols_i=i*(nrcols+1);
      for (k=i;k<nrrows;k++) {
	scale += fabs(u[k_nrcols_i]);
	k_nrcols_i+=nrcols;
      }
      if (scale) {
	k_nrcols_i=i*(nrcols+1);
	for (k=i;k<nrrows;k++) {
	  u[k_nrcols_i] /= scale;
	  s += u[k_nrcols_i]*u[k_nrcols_i];
	  k_nrcols_i+=nrcols;
	}
	f=u[i*nrcols+i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	u[i*nrcols+i]=f-g;
	if (i != nrcols-1) {
	  register int j_i=1; // j_i==j-i;
	  for (j=l;j<nrcols;j++) {
	    k_nrcols_i=i*(nrcols+1);
	    for (s=0.0,k=i;k<nrrows;k++) {
	      s += u[k_nrcols_i]*u[k_nrcols_i+j_i];
	      k_nrcols_i+=nrcols;
	    }
	    f=s/h;
	    k_nrcols_i=i*(nrcols+1);
	    for (k=i;k<nrrows;k++) {
	      u[k_nrcols_i+j_i] += f*u[k_nrcols_i];
	      k_nrcols_i+=nrcols;
	    }
	    j_i++;
	  }
	}
	k_nrcols_i=i*(nrcols+1);
	for (k=i;k<nrrows;k++) {
	  u[k_nrcols_i] *= scale;
	  k_nrcols_i+=nrcols;
	}
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if ((i < nrrows) && (l != nrcols)) {
      register int i_nrcols=i*nrcols;
      for (k=l;k<nrcols;k++) scale += fabs(u[i_nrcols+k]);
      if (scale) {
	for (k=l;k<nrcols;k++) {
	  u[i_nrcols+k] /= scale;
	  s += u[i_nrcols+k]*u[i_nrcols+k];
	}
	f=u[i_nrcols+l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	u[i_nrcols+l]=f-g;
	for (k=l;k<nrcols;k++) rv1[k]=u[i_nrcols+k]/h;
	if (l != nrrows) {
	  register int j_nrcols=l*nrcols;
	  for (j=l;j<nrrows;j++) {
	    for (s=0.0,k=l;k<nrcols;k++) s += u[j_nrcols+k]*u[i_nrcols+k];
	    for (k=l;k<nrcols;k++) u[j_nrcols+k] += s*rv1[k];
	    j_nrcols+=nrcols;
	  }
	}
	for (k=l;k<nrcols;k++) u[i_nrcols+k] *= scale;
      }
    }
    anorm=max(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  register int i_nrcols=nrcols*nrcols;
  for (i=nrcols-1;i>=0;i--) {
    i_nrcols-=nrcols;
    if (l < nrcols) {
      if (g) {
	register int j_nrcols_i=l*nrcols+i;
	register int i_nrcols_l=i_nrcols+l;
	register int j_l=0; // j_l==j-l;
	for (j=l;j<nrcols;j++) {
	  v[j_nrcols_i]=(u[i_nrcols_l+j_l]/u[i_nrcols_l])/g;
	  j_nrcols_i+=nrcols;
	  j_l++;
	}
	register int j_i=1;;
	for (j=l;j<nrcols;j++) {
	  register int k_nrcols_j=l*nrcols+j;
	  for (s=0.0,k=l;k<nrcols;k++) {
	    s += u[i_nrcols+k]*v[k_nrcols_j];
	    k_nrcols_j+=nrcols;
	  }
	  register int k_nrcols_i=l*nrcols+i;
	  for (k=l;k<nrcols;k++) {
	    v[k_nrcols_i+j_i] += s*v[k_nrcols_i];
	    k_nrcols_i+=nrcols;
	  }
	  j_i++;
	}
      }
      register int j_nrcols_i=l*nrcols+i;
      for (j=l;j<nrcols;j++) {
	v[i_nrcols+j]=v[j_nrcols_i]=0.0;
	j_nrcols_i+=nrcols;
      }
    }
    v[i_nrcols+i]=1.0;
    g=rv1[i];
    l=i;
  }
  i_nrcols=nrcols*nrcols;
  for (i=nrcols-1;i>=0;i--) {
    i_nrcols-=nrcols;
    l=i+1;
    g=w[i];
    if (l < nrcols)
      for (j=l;j<nrcols;j++) u[i_nrcols+j]=0.0;
    if (g) {
      g=1.0/g;
      if (l != nrcols) {
	register int j_i=1; // j_i==j-i;
	for (j=l;j<nrcols;j++) {
	  register int k_nrcols_i=l*nrcols+i;
	  for (s=0.0,k=l;k<nrrows;k++) {
	    s += u[k_nrcols_i]*u[k_nrcols_i+j_i];
	    k_nrcols_i+=nrcols;
	  }
	  f=(s/u[i_nrcols+i])*g;
	  k_nrcols_i=i_nrcols+i;
	  for (k=i;k<nrrows;k++) {
	    u[k_nrcols_i+j_i] += f*u[k_nrcols_i];
	    k_nrcols_i+=nrcols;
	  }
	  j_i++;
	}
      }
      register int j_nrcols_i=i_nrcols+i;
      for (j=i;j<nrrows;j++) {
	u[j_nrcols_i] *= g;
	j_nrcols_i+=nrcols;
      }
    }
    else {
      register int j_nrcols_i=i_nrcols+i;
      for (j=i;j<nrrows;j++) {
	u[j_nrcols_i]=0.0;
	j_nrcols_i+=nrcols;
      }
    }
    ++u[i_nrcols+i];
  }
  for (k=nrcols-1;k>=0;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>0;l--) {
	nm=l-1;
	if (fabs(rv1[l])+anorm == anorm) {
	  flag=0;
	  break;
	}
	if (fabs(w[nm])+anorm == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  if (fabs(f)+anorm != anorm) {
	    g=w[i];
	    h=PYTHAG(f,g);
	    w[i]=h;
	    h=1.0/h;
	    c=g*h;
	    s=(-f*h);
	    register int j_nrcols=0;
	    for (j=0;j<nrrows;j++) {
	      y=u[j_nrcols+nm];
	      z=u[j_nrcols+i];
	      u[j_nrcols+nm]=y*c+z*s;
	      u[j_nrcols+i]=z*c-y*s;
	      j_nrcols+=nrcols;
	    }
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  register int j_nrcols_k=k;
	  for (j=0;j<nrcols;j++) {
	    v[j_nrcols_k]=(-v[j_nrcols_k]);
	    j_nrcols_k+=nrcols;
	  }
	}
	break;
      }
      if (its == 30) {
	DL_dsystem->get_companion()->Msg("Warning: DL_largematrix::svdcmp: no convergence in 30 iterations\n");
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=PYTHAG(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=PYTHAG(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y=y*c;
	register int jj_nrcols=0;
	for (jj=0;jj<nrcols;jj++) {
	  x=v[jj_nrcols+j];
	  z=v[jj_nrcols+i];
	  v[jj_nrcols+j]=x*c+z*s;
	  v[jj_nrcols+i]=z*c-x*s;
	  jj_nrcols+=nrcols;
	}
	z=PYTHAG(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=(c*g)+(s*y);
	x=(c*y)-(s*g);
	jj_nrcols=0;
	for (jj=0;jj<nrrows;jj++) {
	  y=u[jj_nrcols+j];
	  z=u[jj_nrcols+i];
	  u[jj_nrcols+j]=y*c+z*s;
	  u[jj_nrcols+i]=z*c-y*s;
	  jj_nrcols+=nrcols;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  delete[] rv1;
  // post process w:
  DL_Scalar biggest=0.0;
  for (i=0;i<nrelem;i+=nrcols+1) if (fabs(a[i])>biggest) biggest=fabs(a[i]);
  // biggest is now the largest magnitude on A's diagonal
#define TINY 1.0e-10
  biggest*=TINY;
#undef TINY
  for (jj=0,i=0;i<nrcols;i++) if (fabs(w[i])<biggest) { w[i]=0.0; jj++; }
if (jj>0) {
  DL_dsystem->get_companion()->Msg("svdcomp: %d singularities\n", jj );
}
  return jj;
}
#undef SIGN
#undef PYTHAG

void DL_largematrix::svbksb(DL_largevector *x, DL_largevector *b){
// solves x from Ax=b using the sv decomposition stored in this matrix
// pre: rep==svdmpd
  int i,j;
  DL_Scalar s, *tmp;

  tmp=new DL_Scalar[nrcols];
  for (j=0;j<nrcols;j++) {
    s=0.0;
    if (w[j]) {
      for (i=0;i<nrcols;i++) s+=u[i*nrcols+j]*b->get(i);
      s/=w[j];
    }
    tmp[j]=s;
  }
  for (j=0;j<nrcols;j++) {
    s=0.0;
    for (i=0;i<nrcols;i++) s+=v[j*nrcols+i]*tmp[i];
    x->set(j,s);
  }
  delete[] tmp;
}

void DL_largematrix::show(){
  int r,c;
  char s[80],t[160]="";
  for (r=0;r<nrrows;r++) {
    for (c=0;c<nrcols;c++) {
      sprintf(s, " %f", a[nrcols*r+c]);
      strcat(t,s);
    }
    strcat(t,"\n");
  }
  DL_dsystem->get_companion()->Msg("%s",t);
}

void DL_largematrix::show_all(){
  int r,c;
  char s[80],t[80];
  DL_dsystem->get_companion()->Msg("matrix:\n");
  show();
  DL_dsystem->get_companion()->Msg("solve method: %s",
				(sm==lud_bcksub ? "LU Decomposition\n" :
				(sm==conjug_grad ? "Conjugate Gradient\n" :
				                    "Singular Value Decomposition\n" ))
			       );
  switch (rep) {
  case full: DL_dsystem->get_companion()->Msg("rep=full\n"); break;
  case riss: DL_dsystem->get_companion()->Msg("rep=riss\n"); break;
  case lud:
  case ludb:
    sprintf(t, "bandw: %d\nlu:\n", bandw );
    for (r=0;r<nrrows;r++) {
      for (c=0;c<nrcols;c++) {
	sprintf(s, " %f", lu[nrcols*r+c]);
	strcat(t,s);
      }
      strcat(t,"\n");
    }
    DL_dsystem->get_companion()->Msg("%s",t);
    break;
  case svdcmpd: DL_dsystem->get_companion()->Msg("rep=svdcmpd\n"); break;
  }
  if (nonzero) {
    sprintf(t,"nonzero:\n");
    for (r=0;r<nrrows;r++) {
      for (c=0;c<nrcols;c++) {
	sprintf(s, " %s", (nonzero[nrcols*r+c]?"X":"-") );
	strcat(t,s);
      }
      strcat(t,"\n");
    }
    DL_dsystem->get_companion()->Msg("%s",t);
  }    
}

void DL_largematrix::set_min_solve_method(solve_method _sm){
  min_sm=_sm;
  set_solve_method(min_sm);
}

void DL_largematrix::set_solve_method(solve_method _sm){
  if (nrcols!=nrrows) {
    DL_dsystem->get_companion()->Msg("Warning: DL_largematrix::set_solve_method():\n Can only solve square systems, and this matrix is not square!!!\n");
    return;
  }
  if (_sm<min_sm) _sm=min_sm;
  if (_sm==sm) return;
#ifdef DEBUG
  DL_dsystem->get_companion()->Msg("Switching from %s solving to %s solving\n",
                                (sm==lud_bcksub ? "LU Decomposition" :
                                (sm==conjug_grad ? "Conjugate Gradient" :
                                                    "Singular Value Decomposition")),
                                (_sm==lud_bcksub ? "LU Decomposition" :
                                (_sm==conjug_grad ? "Conjugate Gradient" :
                                                     "Singular Value Decomposition"))
                               );
#endif
#undef DEBUG
  sm=_sm;
  analyse_structure();
}

void DL_largematrix::analyse_structure(){
  if (nrcols!=nrrows) {
    DL_dsystem->get_companion()->Msg("Warning: DL_largematrix::analyse_structure():\n Can only solve square systems, and this matrix is not square!!!\n");
    return;
  }
  switch (sm) {
  case lud_bcksub:
    get_bandwidth();
    return;
  case conjug_grad:
    if (2*get_nrnonzero()<nrelem) // use sparse matrix representation
      full2riss();
    return;
  case svd:
    return;
  }
}

boolean DL_largematrix::prep_for_solve(){
// returns if there were any singularities
  switch (sm) {
  case lud_bcksub:
    if ((2*bandw>nrrows?ludcmp():ludcmpbw())>=0) {
      // ((near) singular value detected)
      set_solve_method(conjug_grad);
      return prep_for_solve();
    }
    return FALSE;
  case conjug_grad:
    if (2*nrnonzero<nrelem)
      // use previously calculated sparse matrix representation
      rep=riss;
    return FALSE;
  case svd:
    return (svdcmp()!=0);
  }
  return FALSE;
}

boolean DL_largematrix::solve(DL_largevector *x, DL_largevector *b){
// returns if the solve_method has changed
  switch (rep) {
  case full:
  case riss:
    if (conjug_gradient(x,b)) {
      // solution was diverging: so we have a singular matrix and
      // we have to use SVD
      reptofull();
      set_solve_method(svd);
      prep_for_solve();
      solve(x,b);
      return TRUE;
    }
    return FALSE;
  case lud: lubksb(x,b); return FALSE;
  case ludb: lubksbbw(x,b); return FALSE;
  case svdcmpd: svbksb(x,b); return FALSE;
  }
  return FALSE;
}
