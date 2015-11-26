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
// filename: dyna.h
// description: geometries with dynamic motion behaviour
// author: Bart Barenbrug   March 1996
//

#ifndef DL_DYNAH
#define DL_DYNAH

#include "boolean.h"
#include "geo.h"
#include "NaN.h"

// Class Mpair is internal to DL
// Elements of class Mpair are used in the list forces of each dyna
// which administrates pairs of forces/application-points (in
// local coordinates) so torques can be calculated later based on
// the current orientation (so we can have torques that do not rotate
// along with the dyna).

class DL_Mpair : public DL_ListElem {
  public:
    DL_vector force;
    DL_vector rho;
};

// ************* //
// class DL_dyna //
// ************* //

class DL_dyna : public DL_geo {
  protected:
    int     oneD;       // 0 if the object is 2D or 3D; otherwise: the index
                        // of the basevector which is the axis of the dyna.
    DL_vector  J;       // diagonal of diagonalised inertia tensor
    DL_vector  Jinv;    // diagonal of inverse of diagonalised inertia tensor
    DL_vector  F;       // total amount of force applied to this dyna
                        // during this frame
    DL_vector  M;       // total amount central torque applied to this dyna
                        // during this frame (through applytorque)
    DL_supvec mstateimp;    // motion state at cuurent time, plus impulses
    boolean   Fuptodate;    // validity of nextstate wrt F
    boolean   Muptodate;    // validity of nextstate wrt M
    DL_vector FSave;        // F saved for during testing
    DL_vector MSave;        // M saved for during testing
    DL_supvec mstateSave;   // nextmstate saved for during testing
    DL_supvec mstateimpSave;// mstateimp saved for during testing
    boolean   FuptodateSave;// Fuptodate saved for during testing
    boolean   MuptodateSave;// Muptodate saved for during testing
    boolean   testing;      // are we testing?
    DL_Scalar velodamping;  // velocity damping (per time unit)
    DL_Scalar totalmass;    // total mass
    DL_Scalar totalmass_inv;// 1.0/totalmass
    DL_vector Fexternal; // total of external central forces

    DL_List   forces; // reaction force/application point pairs for calculating M
    DL_List   forcesSave;

    // matrix caches for analytical inverse dynamics support:
    // some are not full matrices ((anti)symmetrical), so we
    // only store the relevant elements
    // cache1 actually stores (dAi/ddw)/(h*h) and h*h*(ddw/dM)
    boolean matrixcache1empty;
    DL_Scalar dA0ddw01,dA1ddw01,dA2ddw01;
    DL_Scalar dA0ddw02,dA1ddw02,dA2ddw02;
    DL_Scalar dA0ddw12,dA1ddw12,dA2ddw12;
    DL_Scalar ddwdM00, ddwdM01, ddwdM02, ddwdM11, ddwdM12, ddwdM22;
    boolean matrixcache2empty;
    DL_Scalar m00x,m01x,m02x,  m10x,m11x,m12x,  m20x,m21x,m22x,
              m00y,m01y,m02y,  m10y,m11y,m12y,  m20y,m21y,m22y,
              m00z,m01z,m02z,  m10z,m11z,m12z,  m20z,m21z,m22z;
    
    inline void calc_M(DL_matrix*,DL_vector*);
           // calculate total torque as the sum
           // of M and all torques listed in forces

    inline  void integrate();           // integrate the motion state
    void    update_cache1(void);        // makes sure cache1 is filled
    void    update_cache2(void);        // makes sure cache2 is filled

  public:
    void    assign(DL_dyna*,void*);     // assigment

    virtual void set_position(DL_point*);     // set the position of the dyna
    virtual void set_velocity(DL_vector*);    // set the velocity of the dyna
    virtual void set_orientation(DL_matrix*); // set the orientation of the dyna
    virtual void set_angvelocity(DL_vector*); // set the angular velocity of the dyna
    
    void  set_inertiatensor(DL_Scalar,DL_Scalar,DL_Scalar);
    DL_Scalar get_inertiamoment(int);			  
    void  set_mass(DL_Scalar);            // set total mass
    DL_Scalar get_mass();                 // get total mass
			  
    void  set_velodamping(DL_Scalar vd);       // set velocity damping
    DL_Scalar get_velodamping();               // get velocity damping

    DL_Scalar   kinenergy(void);        // get kinetic energy of the dyna
    DL_Scalar   potenergy(void);        // get potential energy of the dyna
                                        // (wrt to external forces);
    DL_Scalar   totenergy(void);        // get total (kinetic+potential) energy
                                        // of the dyna
					
    virtual void move(DL_point*,DL_matrix*);  // give the geo the new position
                                              // and orientation given by the
			                      // parameters, and update the
				              // (angular) velocity accordingly

    void applyforce(DL_point*, DL_geo*, DL_vector*); // apply a force 
    void applycenterforce(DL_vector*);            // apply a central force 
    void applytorque(DL_vector*);                 // apply a torque
    void applyimpulse(DL_point*, DL_geo*, DL_vector*); // apply an impulse

//////////// For internal use only: /////////////

    void    init();              // initialise attributes
    virtual boolean is_dyna(void) { return TRUE; }

    void prepare_for_next_frame(void);  // shift state and prepare for the next frame;
    void new_frame(void); // start of new frame (after user-code
                          // has been handled)
    
    int     get_oneD();                 // returns oneD

    virtual DL_point* get_next_position(void);     // get the geo's next position
    virtual DL_vector* get_next_velocity(void);    // get the geo's next velocity
    virtual DL_matrix* get_next_orientation(void); // get the geo's next orientation
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

    void imp_toworld(DL_point*, DL_point*); // get worldcoordinates of 1st point (given
                                            // in local coordinates) based on mstateimp
    void imp_toworld(DL_vector*,DL_vector*);
                                            // get worldcoordinates of 1st vector (given
                                            // in local coordinates) based on mstateimp
    void get_impvelocity(DL_point*, DL_vector*);
                                            // get velocity of the point (given in local
                                            // coordinates) based on mstateimp
    void get_impvelocity(DL_vector*,DL_vector*);
                                            // get velocity of the vector (given in local
                                            // coordinates) based on mstateimp
    DL_Scalar   impkinenergy(void);         // get kinetic energy of the dyna after impulses
    DL_Scalar   imppotenergy(void);         // get potential energy of the dyna after impulses
                                        // (wrt to external forces)
    DL_Scalar   imptotenergy(void);         // get total (kinetic+potential) energy
                                        // of the dyna after impulses
    DL_Scalar   newkinenergy(void);         // get new kinetic energy of the dyna
    DL_Scalar   newpotenergy(void);         // get new potential energy of the dyna
                                        // (wrt to external forces)
    DL_Scalar   newtotenergy(void);         // get new total (kinetic+potential) energy
                                        // of the dyna					

    inline void ode(DL_supvec*,DL_supvec*, DL_Scalar);   // motion o.d.e.

    // methods to support empirical determination of dc/dF (for inverse dynamics):
    void begintest();                   // start testing: save F,M,A,nextmstate and uptodate
    void endtest();                     // end testing: restore F,M,A,nextmstate and uptodate

    // methods to support analytical determination of dc/dF (for inverse dynamics):
    inline void dpdfq(DL_point*, DL_point*, DL_matrix*);  // corresponds to applyforce/new_toworld(point)
    inline void dpdF(DL_point*, DL_matrix*);              // corresponds to applycenterforce/new_toworld(point)
    inline void dpdM(DL_point*, DL_matrix*);              // corresponds to applytorque/new_toworld(point)
    inline void dpdi(DL_point*, DL_point*, DL_matrix*);   // corresponds to applyimpulse/new_toworld(point)

    inline void dvdfq(DL_vector*, DL_point*, DL_matrix*); // corresponds to applyforce/new_toworld(vector)
    // zero derivative...                          // corresponds to applycenterforce/new_toworld(vector)
    inline void dvdM(DL_vector*, DL_matrix*);             // corresponds to applytorque/new_toworld(vector)
    inline void dvdi(DL_vector*, DL_point*, DL_matrix*);  // corresponds to applyimpulse/new_toworld(vector)

    inline void ddpdfq(DL_point*, DL_point*, DL_matrix*); // corresponds to applyforce/get_newvelocity(point)
    inline void ddpdF(DL_point*, DL_matrix*);             // corresponds to applycenterforce/get_newvelocity(point)
    inline void ddpdM(DL_point*, DL_matrix*);             // corresponds to applytorque/get_newvelocity(point)
    inline void ddpdi(DL_point*, DL_point*, DL_matrix*);  // corresponds to applyimpulse/get_newvelocity(point)

    inline void ddvdfq(DL_vector*, DL_point*, DL_matrix*);// corresponds to applyforce/get_newvelocity(vector)
    // zero derivative...                          // corresponds to applycenterforce/get_newvelocity(vector)
    inline void ddvdM(DL_vector*, DL_matrix*);            // corresponds to applytorque/get_newvelocity(vector)
    inline void ddvdi(DL_vector*, DL_point*, DL_matrix*); // corresponds to applyimpulse/get_newvelocity(vector)

    // zero derivative...                          // corresponds to applyforce/imp_toworld(point)
    // zero derivative...                          // corresponds to applycenterforce/imp_toworld(point)
    // zero derivative...                          // corresponds to applytorque/imp_toworld(point)
    // zero derivative...                          // corresponds to applyimpulse/imp_toworld(point)

    // zero derivative...                          // corresponds to applyforce/imp_toworld(vector)
    // zero derivative...                          // corresponds to applycenterforce/imp_toworld(vector)
    // zero derivative...                          // corresponds to applytorque/imp_toworld(vector)
    // zero derivative...                          // corresponds to applyimpulse/imp_toworld(vector)

    // zero derivative...                          // corresponds to applyforce/get_impvelocity(point)
    // zero derivative...                          // corresponds to applycenterforce/get_impvelocity(point)
    // zero derivative...                          // corresponds to applytorque/get_impvelocity(point)
    inline void ddptdi(DL_point*, DL_point*, DL_matrix*); // corresponds to applyimpulse/get_impvelocity(point)

    // zero derivative...                          // corresponds to applyforce/get_impvelocity(vector)
    // zero derivative...                          // corresponds to applycenterforce/get_impvelocity(vector)
    // zero derivative...                          // corresponds to applytorque/get_impvelocity(vector)
    inline void ddvtdi(DL_vector*, DL_point*, DL_matrix*);// corresponds to applyimpulse/get_impvelocity(vector)
    
    void reintegrate() { Fuptodate=Muptodate=FALSE; };
         // reintegrate (called by DL_dsystem after integrator is changed)
    DL_Scalar torquefactor(); // returns the factor that a constraint can use to
                          // scale torques to get them like forces in magnitude

               DL_dyna(void*);       // constructor
	       ~DL_dyna();           // destructor
};

#include "dyna_system.h"

inline void DL_dyna::init() {
  testing=FALSE;
  oneD=0;
  J.init(1,1,1);
  Jinv.init(1,1,1);
  F.init(0,0,0);
  M.init(0,0,0);
  matrixcache1empty=matrixcache2empty=TRUE;
  
  Fuptodate=Muptodate=FALSE;

  velodamping=1.0;
  totalmass=totalmass_inv=0.0;

  Fexternal.init(0,0,0);
}

inline DL_dyna::DL_dyna(void *comp):DL_geo(comp) {
  init();
  if (DL_dsystem) {
    DL_dsystem->register_dyna(this);
    DL_dsystem->get_companion()->get_first_geo_info(this);
  }
  else
    fprintf(stderr,"Severe warning: there is no dyna system to manage the dyna's!\n"); // can't use Msg here!!
}

inline DL_dyna::~DL_dyna() {
   if (DL_dsystem)  DL_dsystem->remove_dyna(this);
   forces.delete_all();
   forcesSave.delete_all();
}

inline void DL_dyna::set_position(DL_point *p) {
  DL_geo::set_position(p);
  mstateimp.z.assign(&(mstate.z));
}

inline void DL_dyna::set_velocity(DL_vector *v){
  DL_geo::set_velocity(v);
  mstateimp.v.assign(&(mstate.v));
}

inline void DL_dyna::set_orientation(DL_matrix *m) {
  DL_geo::set_orientation(m);
  mstate.A2q();
  mstate.q2A(this);  
  mstateimp.A.assign(&(mstate.A));
  mstateimp.q.assign(&(mstate.q));  
}

inline void DL_dyna::set_angvelocity(DL_vector *w){
  DL_geo::set_angvelocity(w);
  mstateimp.w.assign(&(mstate.w));
}

inline void DL_dyna::set_velodamping(DL_Scalar vd){
  velodamping=vd;
}

inline DL_Scalar DL_dyna::get_velodamping() {
  return velodamping;
}

inline DL_point* DL_dyna::get_next_position(void) {
  integrate();
  return &(nextmstate.z);
}

inline DL_vector* DL_dyna::get_next_velocity(void) {
  integrate();
  return &(nextmstate.v);
}

inline DL_matrix* DL_dyna::get_next_orientation(void) {
  integrate();
  return &(nextmstate.A);
}

inline DL_vector* DL_dyna::get_next_angvelocity(void) {
  integrate();
  return &(nextmstate.w);
}

inline void DL_dyna::new_toworld(DL_point *p, DL_point *pr) {
  integrate();
  DL_geo::new_toworld(p,pr);
}

inline void DL_dyna::new_toworld(DL_vector *v, DL_vector *vr) {
  integrate();
  DL_geo::new_toworld(v,vr);
}

inline void DL_dyna::new_tolocal(DL_point *p, DL_geo *g, DL_point *pr) {
  if (g!=this) integrate();
  DL_geo::new_tolocal(p,g,pr);
}

inline void DL_dyna::new_tolocal(DL_vector *v, DL_geo *g, DL_vector *vr) {
  if (g!=this) integrate();
  DL_geo::new_tolocal(v,g,vr);
}

inline void DL_dyna::get_newvelocity(DL_point *p, DL_vector *v) {
// p in local coordinates; v in world coordinates
  integrate();
  DL_geo::get_newvelocity(p,v);
}

inline void DL_dyna::get_newvelocity(DL_vector *v, DL_vector *vr) {
// v in local coordinates; vr in world coordinates
  integrate();
  DL_geo::get_newvelocity(v,vr);
}

inline void DL_dyna::imp_toworld(DL_point *p, DL_point *pr) {
  DL_vector vtmp;
  integrate();
  mstateimp.A.times(p,pr);
  mstateimp.z.tovector(&vtmp);
  pr->plusis(&vtmp);
}

inline void DL_dyna::imp_toworld(DL_vector *v, DL_vector *vr) {
  integrate();
  mstateimp.A.times(v,vr);
}

inline void DL_dyna::get_impvelocity(DL_point *p, DL_vector *v) {
// p in local coordinates; v in world coordinates
  DL_point pw;
  DL_vector vw;
  mstateimp.A.times(p,&pw);
  pw.tovector(&vw);
  vw.crossprod(&(mstateimp.w),v);
  v->plusis(&(mstateimp.v));
}

inline void DL_dyna::get_impvelocity(DL_vector *v, DL_vector *vr) {
// v in local coordinates; vr in world coordinates
  DL_vector vw;
  mstateimp.A.times(v,&vw);
  vw.crossprod(&(mstateimp.w),vr);
}

inline int DL_dyna::get_oneD() {
  return oneD;
}

inline void DL_dyna::set_inertiatensor(DL_Scalar j0 ,DL_Scalar j1, DL_Scalar j2) {
  J.init(j0,j1,j2);

  // take care of 1D objects
  {
    oneD=0;
    DL_Scalar temp=J.norm()/1000.0;
    if (fabs(J.x)<temp) oneD=1;
    if (fabs(J.y)<temp) oneD=2;
    if (fabs(J.z)<temp) oneD=3;
  }
  
  // and invert J
  Jinv.x=(oneD!=1 ? 1.0/J.x : 0.0);
  Jinv.y=(oneD!=2 ? 1.0/J.y : 0.0);
  Jinv.z=(oneD!=3 ? 1.0/J.z : 0.0);

  Muptodate=FALSE;
}

inline DL_Scalar DL_dyna::get_inertiamoment(int i) {
  switch (i) {
  case 0: return J.x;
  case 1: return J.y;
  case 2: return J.z;
  }
  return 0;
}

inline void DL_dyna::set_mass(DL_Scalar tm) {
  totalmass=tm;
  if (tm!=0) totalmass_inv=1.0/tm;
  Fuptodate=FALSE;
}
			  
inline DL_Scalar DL_dyna::get_mass() {
  return totalmass;
}

inline void DL_dyna::integrate() {
  if (!Fuptodate) {
    // first do the positional integration:
    // this can be done analytically:
    // z_{t+h}=z_t+v_t*h+0.5*(F/m)*h^2
    // v_{t+h}=v_t+(F/m)*h
    DL_vector vt;
    DL_Scalar h=DL_dsystem->get_integrator()->stepsize();
    F.times(h*totalmass_inv,&(nextmstate.v));
    nextmstate.v.times(0.5,&vt);
    vt.plusis(&(mstateimp.v));
    vt.timesis(h);
    mstateimp.z.plus(&vt,&(nextmstate.z));
    nextmstate.v.plusis(&(mstateimp.v));
    Fuptodate=TRUE;
  }
  if (!Muptodate) {
    // then do the orientational integration
    // using the motion integrator provided by DL_dsystem
    DL_dsystem->get_integrator()->integrate(&mstateimp,this,&nextmstate);
    Muptodate=TRUE;
  }
}

inline void DL_dyna::ode(DL_supvec* y,DL_supvec* dy, DL_Scalar df) {
// this method implements the motion differential equations for orientation
// dw=f(w,A,M) and dq=f(q,w+df*dw)
// df is a discretisation factor used to compensate for a bit of the linearity
// of the euler integrators

  DL_vector v0, v1, koppel; // auxilary variables

  y->A.diag_transpose_vec(&J,&(y->w),&v0); // v0:=J A^T w  
  y->A.times(&v0,&v1);                     // v1:=A J A^T w
  v1.crossprod(&(y->w),&v0);               // v0:=w~ A J A^T w
  calc_M(&(y->A),&koppel);
  koppel.minus(&v0,&v1);	           // v1:=M - w~ A J A^T w
  y->A.diag_transpose_vec(&Jinv,&v1,&v0);  // v0:=Jinv A^T (M - w~ A J A^T w)
  y->A.times(&v0,&(dy->w));                // dy->w:=A Jinv A^T ( M - w~ A J A^T w )
  
  if (df==0.0)
    y->q.times(&(y->w),&(dy->q));  // dy->q:=(y->w)#(y->q)
  else {
    dy->w.times(df,&v0);
    v0.plusis(&(y->w));
    y->q.times(&v0,&(dy->q));      // dy->q:=(y->w + df*dy->w)#(y->q)
  }

}

inline void DL_dyna::calc_M(DL_matrix *A, DL_vector *totkoppel){
  DL_vector t;
  DL_vector tp;
  DL_vector *l;
  DL_Mpair *elem;
  if (oneD!=0) {
    if (oneD==1) l=&(A->c0);
    if (oneD==2) l=&(A->c1);
    if (oneD==3) l=&(A->c2);
  }
  totkoppel->assign(&M);
  elem=(DL_Mpair*)forces.getfirst();
  while (elem) {
    A->times(&(elem->rho),&tp);
    elem->force.crossprod(&tp,&t);
    if (oneD!=0) { // no torque in the direction of l;
      l->times(l->inprod(&t),&tp);
      t.minus(&tp,&t);
    }
    totkoppel->plusis(&t);
    elem=(DL_Mpair*)forces.getnext(elem);
  }
}

inline DL_Scalar DL_dyna::torquefactor() {
  return ::sqrt(J.norm())/totalmass;
}

inline void DL_dyna::update_cache1() {
  // first see if we have to fill the cache or if it is already filled:
  if (matrixcache1empty) {
     // calculate dA0ddw, dA1ddw, dA2ddw and ddwdM
     // we do not fill the diagonals of dAiddw since they are zero
     // anyway and they are antisymmetric too, so we only store the
     // lower halfs. We also shift the factor h*h present in all
     // dAiddw to matrix ddwdM

     // first calculate some of the products:
     DL_Scalar *q=(DL_Scalar*)&(mstate.q.c); // tricky but fast
     DL_Scalar
       q00=q[0]*q[0],
       q01=q[0]*q[1], q11=q[1]*q[1],
       q02=q[0]*q[2], q12=q[1]*q[2], q22=q[2]*q[2],
       q03=q[0]*q[3], q13=q[1]*q[3], q23=q[2]*q[3], q33=q[3]*q[3];

     // use them for the dA_i/ddw:
	
     dA0ddw01=q02+q13;
     dA0ddw02=q03-q12;
     dA0ddw12=0.5*(q00+q11-q22-q33);

     dA1ddw01=q23-q01;
     dA1ddw02=0.5*(q11-q00+q33-q22);
     dA1ddw12=q12+q03;

     dA2ddw01=0.5*(q00-q11-q22+q33);
     dA2ddw02=-(q01+q23);
     dA2ddw12=q13-q02;
      
     // now calculate ddw/dM which is A*diag(Jinv)*A^T which is positive
     // and symmetric
     // account for factor h*h:
     DL_Scalar hh=DL_dsystem->get_integrator()->stepsize(); hh*=hh;
     DL_Scalar Jx=hh*Jinv.x, Jy=hh*Jinv.y, Jz=hh*Jinv.z;
     ddwdM00=Jx*mstate.A.c0.x*mstate.A.c0.x +
             Jy*mstate.A.c1.x*mstate.A.c1.x +
	     Jz*mstate.A.c2.x*mstate.A.c2.x ;
     ddwdM11=Jx*mstate.A.c0.y*mstate.A.c0.y +
             Jy*mstate.A.c1.y*mstate.A.c1.y +
	     Jz*mstate.A.c2.y*mstate.A.c2.y ;
     ddwdM22=Jx*mstate.A.c0.z*mstate.A.c0.z +
             Jy*mstate.A.c1.z*mstate.A.c1.z +
	     Jz*mstate.A.c2.z*mstate.A.c2.z ;
     ddwdM01=Jx*mstate.A.c0.x*mstate.A.c0.y +
             Jy*mstate.A.c1.x*mstate.A.c1.y +
	     Jz*mstate.A.c2.x*mstate.A.c2.y ;
     ddwdM02=Jx*mstate.A.c0.x*mstate.A.c0.z +
             Jy*mstate.A.c1.x*mstate.A.c1.z +
	     Jz*mstate.A.c2.x*mstate.A.c2.z ;
     ddwdM12=Jx*mstate.A.c0.y*mstate.A.c0.z +
             Jy*mstate.A.c1.y*mstate.A.c1.z +
	     Jz*mstate.A.c2.y*mstate.A.c2.z ;
      
     matrixcache1empty=FALSE;
  }
}

inline void DL_dyna::dpdfq(DL_point *p, DL_point *q, DL_matrix *m) {
// what is the effect on point p if we exert a force to point q
// (both p and q in local coordinates)
  DL_matrix dpdm;
  DL_point rho;
  DL_Scalar hhm=DL_dsystem->get_integrator()->stepsize();
  hhm*=0.5*hhm*totalmass_inv; // 0.5*h*h/m
    
  dpdM(p,&dpdm);
  // we now have dp/dM. we know dM/dfq=((mstate.A)*q)~
  // so we multiply the two and add the bit for dp/dF

  // first calculate (mstate.A)*q :
  mstate.A.times(q,&rho);
  
  // now multiply: we can do this a bit cheaper since we know the
  // diagonal of rho~ only contains zeros. Of course I should not
  // really use the implementations of matrices and vectors like
  // this, but these routines need to be FAST:
  m->c0.x=dpdm.c2.x*rho.y-dpdm.c1.x*rho.z+hhm;
  m->c1.x=dpdm.c0.x*rho.z-dpdm.c2.x*rho.x;
  m->c2.x=dpdm.c1.x*rho.x-dpdm.c0.x*rho.y;
  m->c0.y=dpdm.c2.y*rho.y-dpdm.c1.y*rho.z;
  m->c1.y=dpdm.c0.y*rho.z-dpdm.c2.y*rho.x+hhm;
  m->c2.y=dpdm.c1.y*rho.x-dpdm.c0.y*rho.y;
  m->c0.z=dpdm.c2.z*rho.y-dpdm.c1.z*rho.z;
  m->c1.z=dpdm.c0.z*rho.z-dpdm.c2.z*rho.x;
  m->c2.z=dpdm.c1.z*rho.x-dpdm.c0.z*rho.y+hhm;
}

inline void DL_dyna::dvdfq(DL_vector *v, DL_point *q, DL_matrix *m) {
// what is the effect on vector v if we exert a force to point q
// (both v and q in local coordinates)
// based on v=A*rho; see dpdfq which is based on p=z+A*rho
  DL_matrix dvdm;
  DL_point rho;
  dvdM(v,&dvdm);
  mstate.A.times(q,&rho);
  m->c0.x=dvdm.c2.x*rho.y-dvdm.c1.x*rho.z;
  m->c1.x=dvdm.c0.x*rho.z-dvdm.c2.x*rho.x;
  m->c2.x=dvdm.c1.x*rho.x-dvdm.c0.x*rho.y;
  m->c0.y=dvdm.c2.y*rho.y-dvdm.c1.y*rho.z;
  m->c1.y=dvdm.c0.y*rho.z-dvdm.c2.y*rho.x;
  m->c2.y=dvdm.c1.y*rho.x-dvdm.c0.y*rho.y;
  m->c0.z=dvdm.c2.z*rho.y-dvdm.c1.z*rho.z;
  m->c1.z=dvdm.c0.z*rho.z-dvdm.c2.z*rho.x;
  m->c2.z=dvdm.c1.z*rho.x-dvdm.c0.z*rho.y;
}

inline void DL_dyna::dpdF(DL_point *p, DL_matrix *m) {
// what is the effect on point p if we exert a central force
// on the dyna
  DL_Scalar h=DL_dsystem->get_integrator()->stepsize();
  m->c0.x=m->c1.y=m->c2.z=0.5*h*h*totalmass_inv;
  m->c0.y=m->c0.z=m->c1.x=m->c1.z=m->c2.x=m->c2.y=0.0;
}

inline void DL_dyna::dpdM(DL_point *p, DL_matrix *m) {
// what is the effect on point p if we exert a central torque
// on the dyna

// dp/dM=(Si: i=0,1,2: p_i * (dA_i/ddw)) * (ddw/dM)
// we have the expressions for (dA_i/ddw) and (ddw/dM) is the
// inertia tensor in world coordinates. The (dA_i/ddw) en (ddw/dM)
// matrices can be shared between different calls to this method
// in the same timeframe, so we cache them in variables dA0ddw,
// dA1ddw,dA2ddw and ddwdM (note that is a lot cheaper than caching
// matrices (dA_i/dM) since the three matrix-multiplications that
// are then required initialy probably are more exspensive than the
// one required if (dA_i/ddw) and (ddw/dM) are cached (usually the
// matrices are shared between two calls, sometimes only one, and
// only very rarely three or more.

   // make sure the cache is up to date
   update_cache1();
   
   // now combine the dAiddw matrices with the given p to
   // arrive at dp/ddw. Again this matrix is antisymmetrical
   // with zeros on the diagonal, so we only have to compute
   // three of its elements

   DL_Scalar dpddw01, dpddw02, dpddw12;
   dpddw01=p->x*dA0ddw01 + p->y*dA1ddw01 + p->z*dA2ddw01;
   dpddw02=p->x*dA0ddw02 + p->y*dA1ddw02 + p->z*dA2ddw02;
   dpddw12=p->x*dA0ddw12 + p->y*dA1ddw12 + p->z*dA2ddw12;

   // now multiply dp/ddw with ddw/dM to arrive at the
   // end result dp/dM
   DL_Scalar tmp0=dpddw01*ddwdM01;
   DL_Scalar tmp1=dpddw02*ddwdM02;
   DL_Scalar tmp2=dpddw12*ddwdM12;
   m->c0.x=-(tmp0+tmp1);
    m->c1.x=-(dpddw01*ddwdM11+dpddw02*ddwdM12);
      m->c2.x=-(dpddw01*ddwdM12+dpddw02*ddwdM22);
   m->c0.y=dpddw01*ddwdM00-dpddw12*ddwdM02;
    m->c1.y=tmp0-tmp2;
      m->c2.y=dpddw01*ddwdM02-dpddw12*ddwdM22;
   m->c0.z=dpddw02*ddwdM00+dpddw12*ddwdM01;
    m->c1.z=dpddw02*ddwdM01+dpddw12*ddwdM11;
      m->c2.z=tmp1+tmp2;   
}

inline void DL_dyna::dvdM(DL_vector *v, DL_matrix *m) {
// what is the effect on vector v if we exert a central torque
// on the dyna (see dpdM)

   update_cache1();
   
   DL_Scalar dpddw01, dpddw02, dpddw12;
   dpddw01=v->x*dA0ddw01 + v->y*dA1ddw01 + v->z*dA2ddw01;
   dpddw02=v->x*dA0ddw02 + v->y*dA1ddw02 + v->z*dA2ddw02;
   dpddw12=v->x*dA0ddw12 + v->y*dA1ddw12 + v->z*dA2ddw12;

   DL_Scalar tmp0=dpddw01*ddwdM01;
   DL_Scalar tmp1=dpddw02*ddwdM02;
   DL_Scalar tmp2=dpddw12*ddwdM12;
   m->c0.x=-(tmp0+tmp1);
    m->c1.x=-(dpddw01*ddwdM11+dpddw02*ddwdM12);
      m->c2.x=-(dpddw01*ddwdM12+dpddw02*ddwdM22);
   m->c0.y=dpddw01*ddwdM00-dpddw12*ddwdM02;
    m->c1.y=tmp0-tmp2;
      m->c2.y=dpddw01*ddwdM02-dpddw12*ddwdM22;
   m->c0.z=dpddw02*ddwdM00+dpddw12*ddwdM01;
    m->c1.z=dpddw02*ddwdM01+dpddw12*ddwdM11;
      m->c2.z=tmp1+tmp2;   
}

inline void DL_dyna::update_cache2() {
  if (matrixcache2empty) {

    update_cache1();  // results from cache 1 are required here
    DL_vector Ai;
    DL_Scalar hinv=1.0/DL_dsystem->get_integrator()->stepsize();
        
    mstate.A.c0.times(hinv,&Ai);
    DL_Scalar wzc01=mstate.w.z*dA0ddw01;
    DL_Scalar wyc02=mstate.w.y*dA0ddw02;
    DL_Scalar wxc12=mstate.w.x*dA0ddw12;
    
    m00x= wzc01-wyc02;
    m00y= mstate.w.x*dA0ddw02+Ai.z;
    m00z=-mstate.w.x*dA0ddw01-Ai.y;
    m01x=-mstate.w.y*dA0ddw12-Ai.z;
    m01y= wzc01+wxc12;
    m01z=-mstate.w.y*dA0ddw01+Ai.x;
    m02x=-mstate.w.z*dA0ddw12+Ai.y;
    m02y= mstate.w.z*dA0ddw02-Ai.x;
    m02z= wxc12-wyc02;
    
    mstate.A.c1.times(hinv,&Ai);
    wzc01=mstate.w.z*dA1ddw01;
    wyc02=mstate.w.y*dA1ddw02;
    wxc12=mstate.w.x*dA1ddw12;
    
    m10x= wzc01-wyc02;
    m10y= mstate.w.x*dA1ddw02+Ai.z;
    m10z=-mstate.w.x*dA1ddw01-Ai.y;
    m11x=-mstate.w.y*dA1ddw12-Ai.z;
    m11y= wzc01+wxc12;
    m11z=-mstate.w.y*dA1ddw01+Ai.x;
    m12x=-mstate.w.z*dA1ddw12+Ai.y;
    m12y= mstate.w.z*dA1ddw02-Ai.x;
    m12z= wxc12-wyc02;
    
    mstate.A.c2.times(hinv,&Ai);
    wzc01=mstate.w.z*dA2ddw01;
    wyc02=mstate.w.y*dA2ddw02;
    wxc12=mstate.w.x*dA2ddw12;
    
    m20x= wzc01-wyc02;
    m20y= mstate.w.x*dA2ddw02+Ai.z;
    m20z=-mstate.w.x*dA2ddw01-Ai.y;
    m21x=-mstate.w.y*dA2ddw12-Ai.z;
    m21y= wzc01+wxc12;
    m21z=-mstate.w.y*dA2ddw01+Ai.x;
    m22x=-mstate.w.z*dA2ddw12+Ai.y;
    m22y= mstate.w.z*dA2ddw02-Ai.x;
    m22z= wxc12-wyc02;
    
    matrixcache2empty=FALSE;
  }
}

inline void DL_dyna::ddpdfq(DL_point *p, DL_point *q, DL_matrix *m) {
// what is the effect on the velocity of point p if we exert a force to point q
// (both p and q in local coordinates)
  DL_matrix ddpdm;
  DL_point rho;
  DL_Scalar hm=DL_dsystem->get_integrator()->stepsize()*totalmass_inv; // h/m
    
  ddpdM(p,&ddpdm);
  // we now have ddp/dM. we know dM/dfq=((mstate->A)*q)~
  // so we multiply the two and add the bit for ddp/dF

  // first calculate (mstate->A)*q :
  mstate.A.times(q,&rho);
  
  // now multiply: we can do this a bit cheaper since we know the
  // diagonal of rho~ only contains zero's. Of course I shouldn't
  // really use the implementations of matrices and vectors like
  // this, but these routines need to be FAST:
  m->c0.x=ddpdm.c2.x*rho.y-ddpdm.c1.x*rho.z+hm;
  m->c1.x=ddpdm.c0.x*rho.z-ddpdm.c2.x*rho.x;
  m->c2.x=ddpdm.c1.x*rho.x-ddpdm.c0.x*rho.y;
  m->c0.y=ddpdm.c2.y*rho.y-ddpdm.c1.y*rho.z;
  m->c1.y=ddpdm.c0.y*rho.z-ddpdm.c2.y*rho.x+hm;
  m->c2.y=ddpdm.c1.y*rho.x-ddpdm.c0.y*rho.y;
  m->c0.z=ddpdm.c2.z*rho.y-ddpdm.c1.z*rho.z;
  m->c1.z=ddpdm.c0.z*rho.z-ddpdm.c2.z*rho.x;
  m->c2.z=ddpdm.c1.z*rho.x-ddpdm.c0.z*rho.y+hm;
}

inline void DL_dyna::ddvdfq(DL_vector *v, DL_point *q, DL_matrix *m) {
// what is the effect on the velocity of vector v if we exert a force to point q
// (both p and q in local coordinates)
  DL_matrix ddvdm;
  DL_point rho;
    
  ddvdM(v,&ddvdm);
  mstate.A.times(q,&rho);
  
  m->c0.x=ddvdm.c2.x*rho.y-ddvdm.c1.x*rho.z;
  m->c1.x=ddvdm.c0.x*rho.z-ddvdm.c2.x*rho.x;
  m->c2.x=ddvdm.c1.x*rho.x-ddvdm.c0.x*rho.y;
  m->c0.y=ddvdm.c2.y*rho.y-ddvdm.c1.y*rho.z;
  m->c1.y=ddvdm.c0.y*rho.z-ddvdm.c2.y*rho.x;
  m->c2.y=ddvdm.c1.y*rho.x-ddvdm.c0.y*rho.y;
  m->c0.z=ddvdm.c2.z*rho.y-ddvdm.c1.z*rho.z;
  m->c1.z=ddvdm.c0.z*rho.z-ddvdm.c2.z*rho.x;
  m->c2.z=ddvdm.c1.z*rho.x-ddvdm.c0.z*rho.y;
}

inline void DL_dyna::ddpdF(DL_point *p, DL_matrix *m){
// what is the effect on the velocity of point p if we exert a central force
// on the dyna
  DL_Scalar h=DL_dsystem->get_integrator()->stepsize();
  m->c0.x=m->c1.y=m->c2.z=h*totalmass_inv;
  m->c0.y=m->c0.z=m->c1.x=m->c1.z=m->c2.x=m->c2.y=0.0;
}

inline void DL_dyna::ddpdM(DL_point *p, DL_matrix *m){
// what is the effect on the velocity of point p if we exert a central torque
// on the dyna
   DL_Scalar mtmp0x, mtmp1x, mtmp2x,
         mtmp0y, mtmp1y, mtmp2y,
         mtmp0z, mtmp1z, mtmp2z;
	 
   // make sure the cache is up to date:
   update_cache2();

   // then sum it all up:
   mtmp0x=m00x*p->x+m10x*p->y+m20x*p->z;
   mtmp0y=m00y*p->x+m10y*p->y+m20y*p->z;
   mtmp0z=m00z*p->x+m10z*p->y+m20z*p->z;
   mtmp1x=m01x*p->x+m11x*p->y+m21x*p->z;
   mtmp1y=m01y*p->x+m11y*p->y+m21y*p->z;
   mtmp1z=m01z*p->x+m11z*p->y+m21z*p->z;
   mtmp2x=m02x*p->x+m12x*p->y+m22x*p->z;
   mtmp2y=m02y*p->x+m12y*p->y+m22y*p->z;
   mtmp2z=m02z*p->x+m12z*p->y+m22z*p->z;

   //and multiply with ddw/dM
   m->c0.x=mtmp0x*ddwdM00+mtmp1x*ddwdM01+mtmp2x*ddwdM02;
   m->c0.y=mtmp0y*ddwdM00+mtmp1y*ddwdM01+mtmp2y*ddwdM02;
   m->c0.z=mtmp0z*ddwdM00+mtmp1z*ddwdM01+mtmp2z*ddwdM02;
   m->c1.x=mtmp0x*ddwdM01+mtmp1x*ddwdM11+mtmp2x*ddwdM12;
   m->c1.y=mtmp0y*ddwdM01+mtmp1y*ddwdM11+mtmp2y*ddwdM12;
   m->c1.z=mtmp0z*ddwdM01+mtmp1z*ddwdM11+mtmp2z*ddwdM12;
   m->c2.x=mtmp0x*ddwdM02+mtmp1x*ddwdM12+mtmp2x*ddwdM22;
   m->c2.y=mtmp0y*ddwdM02+mtmp1y*ddwdM12+mtmp2y*ddwdM22;
   m->c2.z=mtmp0z*ddwdM02+mtmp1z*ddwdM12+mtmp2z*ddwdM22;
}

inline void DL_dyna::ddvdM(DL_vector *v, DL_matrix *m){
// what is the effect on the velocity of vector v if we exert a central torque
// on the dyna (see ddpdM)
   DL_Scalar mtmp0x, mtmp1x, mtmp2x,
         mtmp0y, mtmp1y, mtmp2y,
         mtmp0z, mtmp1z, mtmp2z;
	 
   update_cache2();

   mtmp0x=m00x*v->x+m10x*v->y+m20x*v->z;
   mtmp0y=m00y*v->x+m10y*v->y+m20y*v->z;
   mtmp0z=m00z*v->x+m10z*v->y+m20z*v->z;
   mtmp1x=m01x*v->x+m11x*v->y+m21x*v->z;
   mtmp1y=m01y*v->x+m11y*v->y+m21y*v->z;
   mtmp1z=m01z*v->x+m11z*v->y+m21z*v->z;
   mtmp2x=m02x*v->x+m12x*v->y+m22x*v->z;
   mtmp2y=m02y*v->x+m12y*v->y+m22y*v->z;
   mtmp2z=m02z*v->x+m12z*v->y+m22z*v->z;

   //and multiply with ddw/dM
   m->c0.x=mtmp0x*ddwdM00+mtmp1x*ddwdM01+mtmp2x*ddwdM02;
   m->c0.y=mtmp0y*ddwdM00+mtmp1y*ddwdM01+mtmp2y*ddwdM02;
   m->c0.z=mtmp0z*ddwdM00+mtmp1z*ddwdM01+mtmp2z*ddwdM02;
   m->c1.x=mtmp0x*ddwdM01+mtmp1x*ddwdM11+mtmp2x*ddwdM12;
   m->c1.y=mtmp0y*ddwdM01+mtmp1y*ddwdM11+mtmp2y*ddwdM12;
   m->c1.z=mtmp0z*ddwdM01+mtmp1z*ddwdM11+mtmp2z*ddwdM12;
   m->c2.x=mtmp0x*ddwdM02+mtmp1x*ddwdM12+mtmp2x*ddwdM22;
   m->c2.y=mtmp0y*ddwdM02+mtmp1y*ddwdM12+mtmp2y*ddwdM22;
   m->c2.z=mtmp0z*ddwdM02+mtmp1z*ddwdM12+mtmp2z*ddwdM22;
}

inline void DL_dyna::dpdi(DL_point *p, DL_point *q, DL_matrix *m) {
// what is the effect on point p (at t+h) if we exert an impulse to point q
// (both p and q in local coordinates)

   // first calculate ddq/di:
   DL_matrix Atmp,ddqdi;
   Atmp.diagcrosstranspose(&Jinv,q,&(mstate.A));
   mstate.A.times(&Atmp,&ddqdi);

   // make sure the cache is up to date
   update_cache1();
   
   // now combine the dAiddw matrices with the given p to
   // arrive at dp/ddq. Again this matrix is antisymmetrical
   // with zeros on the diagonal, so we only have to compute
   // three of its elements

   DL_Scalar dpddq01, dpddq02, dpddq12, fac=2.0*DL_dsystem->get_integrator()->stepsize();
   dpddq01=fac*(p->x*dA0ddw01 + p->y*dA1ddw01 + p->z*dA2ddw01);
   dpddq02=fac*(p->x*dA0ddw02 + p->y*dA1ddw02 + p->z*dA2ddw02);
   dpddq12=fac*(p->x*dA0ddw12 + p->y*dA1ddw12 + p->z*dA2ddw12);

   // now multiply dp/ddq with ddq/di to arrive at the
   // end result dp/di, accounting for the translational part
   fac=DL_dsystem->get_integrator()->stepsize()*totalmass_inv;
   m->c0.x=dpddq01*ddqdi.c0.y+dpddq02*ddqdi.c0.z+fac;
   m->c1.x=dpddq01*ddqdi.c1.y+dpddq02*ddqdi.c1.z;
   m->c2.x=dpddq01*ddqdi.c2.y+dpddq02*ddqdi.c2.z;
   m->c0.y=dpddq12*ddqdi.c0.z-dpddq01*ddqdi.c0.x;
   m->c1.y=dpddq12*ddqdi.c1.z-dpddq01*ddqdi.c1.x+fac;
   m->c2.y=dpddq12*ddqdi.c2.z-dpddq01*ddqdi.c2.x;
   m->c0.z=-dpddq02*ddqdi.c0.x-dpddq12*ddqdi.c0.y;
   m->c1.z=-dpddq02*ddqdi.c1.x-dpddq12*ddqdi.c1.y;
   m->c2.z=-dpddq02*ddqdi.c2.x-dpddq12*ddqdi.c2.y+fac;

}

inline void DL_dyna::dvdi(DL_vector *v, DL_point *q, DL_matrix *m) {
// what is the effect on vector v (at t+h) if we exert an impulse to point q
// (both v and q in local coordinates)

   // first calculate ddq/di:
   DL_matrix Atmp,ddqdi;
   Atmp.diagcrosstranspose(&Jinv,q,&(mstate.A));
   mstate.A.times(&Atmp,&ddqdi);

   // make sure the cache is up to date
   update_cache1();
   
   // now combine the dAiddw matrices with the given p to
   // arrive at dp/ddq. Again this matrix is antisymmetrical
   // with zeros on the diagonal, so we only have to compute
   // three of its elements

   DL_Scalar dvddq01, dvddq02, dvddq12, fac=2.0*DL_dsystem->get_integrator()->stepsize();
   dvddq01=fac*(v->x*dA0ddw01 + v->y*dA1ddw01 + v->z*dA2ddw01);
   dvddq02=fac*(v->x*dA0ddw02 + v->y*dA1ddw02 + v->z*dA2ddw02);
   dvddq12=fac*(v->x*dA0ddw12 + v->y*dA1ddw12 + v->z*dA2ddw12);

   // now multiply dv/ddq with ddq/di to arrive at the
   // end result dv/di
   m->c0.x=dvddq01*ddqdi.c0.y+dvddq02*ddqdi.c0.z;
   m->c1.x=dvddq01*ddqdi.c1.y+dvddq02*ddqdi.c1.z;
   m->c2.x=dvddq01*ddqdi.c2.y+dvddq02*ddqdi.c2.z;
   m->c0.y=dvddq12*ddqdi.c0.z-dvddq01*ddqdi.c0.x;
   m->c1.y=dvddq12*ddqdi.c1.z-dvddq01*ddqdi.c1.x;
   m->c2.y=dvddq12*ddqdi.c2.z-dvddq01*ddqdi.c2.x;
   m->c0.z=-dvddq02*ddqdi.c0.x-dvddq12*ddqdi.c0.y;
   m->c1.z=-dvddq02*ddqdi.c1.x-dvddq12*ddqdi.c1.y;
   m->c2.z=-dvddq02*ddqdi.c2.x-dvddq12*ddqdi.c2.y;
}

inline void DL_dyna::ddpdi(DL_point *p, DL_point *q, DL_matrix *m) {
// what is the effect on the velocity of point p (at t+h) if we exert an impulse to point q
// (both p and q in local coordinates)
		  
    DL_matrix Atmp,dwdi;
    DL_Scalar S01, S02, S12;
    DL_point rho;
    
    update_cache1();  // results from cache 1 are required here

    // first calculate dw/di:
    Atmp.diagcrosstranspose(&Jinv,q,&(mstate.A));
    mstate.A.times(&Atmp,&dwdi);

    // then calculate ddp/dw in Atmp:
    p->times(2.0*DL_dsystem->get_integrator()->stepsize(),&rho);
    S01=rho.x*dA0ddw01+rho.y*dA1ddw01+rho.z*dA2ddw01;
    S02=rho.x*dA0ddw02+rho.y*dA1ddw02+rho.z*dA2ddw02;
    S12=rho.x*dA0ddw12+rho.y*dA1ddw12+rho.z*dA2ddw12;

    mstate.A.times(p,&rho);
    
    Atmp.c0.x= mstate.w.z*S01-mstate.w.y*S02;
    Atmp.c1.x=-mstate.w.y*S12+rho.z;
    Atmp.c2.x=-mstate.w.z*S12-rho.y;
    Atmp.c0.y=-mstate.w.x*S02-rho.z;
    Atmp.c1.y= mstate.w.z*S01+mstate.w.x*S12;
    Atmp.c2.y=-mstate.w.z*S02+rho.x;
    Atmp.c0.z=-mstate.w.x*S01+rho.y;
    Atmp.c1.z=-mstate.w.y*S01-rho.x;
    Atmp.c2.z= mstate.w.x*S12-mstate.w.y*S02;

   // multiply the two and account for the translational part:
   Atmp.times(&dwdi,m);
   m->c0.x+=totalmass_inv; m->c1.y+=totalmass_inv; m->c2.z+=totalmass_inv;
}

inline void DL_dyna::ddvdi(DL_vector *v, DL_point *q, DL_matrix *m) {
// what is the effect on the velocity of vector v (at t+h) if we exert an impulse to point q
// (both v and q in local coordinates)

    DL_matrix Atmp,dwdi;
    DL_Scalar S01, S02, S12;
    DL_vector rho;
    
    update_cache1();  // results from cache 1 are required here

    // first calculate dw/di:
    Atmp.diagcrosstranspose(&Jinv,q,&(mstate.A));
    mstate.A.times(&Atmp,&dwdi);

    // then calculate ddp/dw in Atmp:
    v->times(2.0*DL_dsystem->get_integrator()->stepsize(),&rho);
    S01=rho.x*dA0ddw01+rho.y*dA1ddw01+rho.z*dA2ddw01;
    S02=rho.x*dA0ddw02+rho.y*dA1ddw02+rho.z*dA2ddw02;
    S12=rho.x*dA0ddw12+rho.y*dA1ddw12+rho.z*dA2ddw12;

    mstate.A.times(v,&rho);
    
    Atmp.c0.x= mstate.w.z*S01-mstate.w.y*S02;
    Atmp.c1.x=-mstate.w.y*S12+rho.z;
    Atmp.c2.x=-mstate.w.z*S12-rho.y;
    Atmp.c0.y=-mstate.w.x*S02-rho.z;
    Atmp.c1.y= mstate.w.z*S01+mstate.w.x*S12;
    Atmp.c2.y=-mstate.w.z*S02+rho.x;
    Atmp.c0.z=-mstate.w.x*S01+rho.y;
    Atmp.c1.z=-mstate.w.y*S01-rho.x;
    Atmp.c2.z= mstate.w.x*S12-mstate.w.y*S02;

   // multiply the two:
   Atmp.times(&dwdi,m);
}

inline void DL_dyna::ddptdi(DL_point *p, DL_point *q, DL_matrix *m) {
// what is the effect on the velocity of point p (at t, after impulses) if we exert an impulse to point q
// (both p and q in local coordinates)
  DL_matrix Atmp;
  m->negcrossdiagcross(p,&Jinv,q);
  m->timestranspose(&(mstate.A),&Atmp);
  mstate.A.times(&Atmp,m);
  m->c0.x+=totalmass_inv;
  m->c1.y+=totalmass_inv;
  m->c2.z+=totalmass_inv;
}

inline void DL_dyna::ddvtdi(DL_vector *v, DL_point *q, DL_matrix *m) {
// what is the effect on the velocity of vector v (at t, after impulses) if we exert an impulse to point q
// (both v and q in local coordinates)
  DL_matrix Atmp;
  DL_point p;
  v->topoint(&p);
  m->negcrossdiagcross(&p,&Jinv,q);
  m->timestranspose(&(mstate.A),&Atmp);
  mstate.A.times(&Atmp,m);
}

inline void DL_dyna::move(DL_point *newpos, DL_matrix *neworient){
  DL_geo::move(newpos,neworient);
  mstate.q.assign(&nextmstate.q);
  nextmstate.A2q();
  nextmstate.q2A(this);  
  mstateimp.assign(&mstate);
}

inline void DL_dyna::new_frame(void) {
  Fexternal.assign(&F);
}

inline void DL_dyna::applycenterforce(DL_vector *f) {
  if ((f->x==0.0)&&(f->y==0.0)&&(f->z==0.0)) return;
  F.plusis(f);
  Fuptodate=FALSE;
}

inline void DL_dyna::applyforce(DL_point *p, DL_geo *g, DL_vector *f) {
  DL_vector m;
  DL_vector t;
  DL_Mpair *forceselem;
  if ((f->x==0.0)&&(f->y==0.0)&&(f->z==0.0)) return;
  
  F.plusis(f);
  Fuptodate=FALSE;  

  if (g==this) p->tovector(&m);
  else {
    DL_point pm;
    to_local(p,g,&pm);
    pm.tovector(&m);
  }

  // add the force to the force list:
  // see if we already have a force to this point;
  forceselem=(DL_Mpair*)forces.getfirst();
  while (forceselem) {
    if (forceselem->rho.equal(&m)) break;
    else forceselem=(DL_Mpair*)forces.getnext(forceselem);
  }
  if (forceselem) {
    forceselem->force.plusis(f);
/* probably won't make things faster:
    if (forceselem->force.norm()==0) {
      forces.remelem(forceselem);
      delete forceselem;
    }
*/
  }
  else {
    forceselem=new DL_Mpair;
    forceselem->force.assign(f);
    forceselem->rho.assign(&m);
    forces.addelem(forceselem);
  }

  Muptodate=FALSE;
}

inline void DL_dyna::applytorque(DL_vector *t) {
  if ((t->x==0.0)&&(t->y==0.0)&&(t->z==0.0)) return;
  if (oneD!=0) {
    DL_vector *l;
    if (oneD==1) l=&(mstate.A.c0);
    if (oneD==2) l=&(mstate.A.c1);
    if (oneD==3) l=&(mstate.A.c2);
    DL_vector tp;
    l->times(l->inprod(t),&tp);
    t->minus(&tp,&tp); // tp is the projection of t onto the plane
                       // perpendicular to the direction of the object
    M.plusis(&tp);
  }
  else M.plusis(t);
  Muptodate=FALSE;
}

inline void DL_dyna::applyimpulse(DL_point *p, DL_geo *g, DL_vector *i) {
  DL_vector delta,tmp;
  DL_vector r; // A(p in lc)
  DL_point pm;

  i->times(totalmass_inv,&delta);
  mstateimp.v.plusis(&delta);

  if (g==this) {
    mstate.A.times(p,&pm);
    pm.tovector(&r);
  }
  else {
    g->to_world(p,&pm);
    pm.minus(&(mstate.z),&r);
  }
  i->crossprod(&r,&delta);
  mstate.A.diag_transpose_vec(&Jinv,&delta,&tmp);
  mstate.A.times(&tmp,&delta);
    
  if (oneD!=0) { // no torque in the direction of the oneDth base vector;
    if (oneD==1) mstate.A.c0.times(mstate.A.c0.inprod(&delta),&tmp);
    if (oneD==2) mstate.A.c1.times(mstate.A.c1.inprod(&delta),&tmp);
    if (oneD==3) mstate.A.c2.times(mstate.A.c2.inprod(&delta),&tmp);
    delta.minusis(&tmp);
  }
  
  mstateimp.w.plusis(&delta);

  Fuptodate=FALSE;  
  Muptodate=FALSE;
}

#endif
