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
// filename     : multibar.cpp
// description	: non-inline methods of class DL_multi_bar
// author	: Bart Barenbrug     April '97
//

#include "dyna_system.h"
#include "multibar.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_multi_bar::DL_multi_bar():DL_constraint() {
  dim=1;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  maxforce=0.0;
  nr=0;
  g=NULL; p=NULL; pw=NULL; dpw=NULL; d=NULL;
  l=1.0;
  rope=FALSE;
}

DL_multi_bar::~DL_multi_bar() {
  delete[] g;
  delete[] p; delete[] pw; delete[] dpw; delete[] d;
  delete[] g_is_dyna;
}

void DL_multi_bar::init(int _maxnr) {
  maxnr=_maxnr;
  g_is_dyna=new boolean[maxnr];
  for (int i=0;i<maxnr;i++) g_is_dyna[i]=FALSE;
  g=new DL_geo*[maxnr];
  p=new DL_point[maxnr];
  pw=new DL_point[maxnr];
  dpw=new DL_vector[maxnr];
  d=new DL_vector[maxnr-1];
    
  F->makezero(); oldF->makezero();

  rp.init(this);
}

void DL_multi_bar::set_length(DL_Scalar _l){
  if (_l>0) l=_l;
  else {
    DL_dsystem->get_companion()->Msg("warning multibar::set_length: length should be greater than zero! Old length kept\n");
  }
}

void DL_multi_bar::addpair(DL_geo *_g, DL_point *_p){
  if (nr==maxnr) {
    DL_dsystem->get_companion()->Msg("Warning: DL_multi_bar::addpair: Too many pairs. Pair not added\n");
    return;
  }

  g[nr]=_g;
  p[nr]=_p;
  
  if (_g) g_is_dyna[nr]=_g->is_dyna();
  if (_g) {
    _g->new_toworld(_p,&(pw[nr]));
    _g->get_newvelocity(_p,&(dpw[nr]));
  }
  else {
    pw[nr].assign(_p);
    dpw[nr].init(0,0,0);
  }
  
  nr++;
  if (nr>=2) {
    if (nr==2) DL_constraint::init();
    set_length(actual_length());
  }
}
  
void DL_multi_bar::reactionforce(int i, DL_vector *f) {
  if ((i<0) || (i+1>=nr)) {
    f->init(0,0,0);
    return;
  }
  
  d[i].times(F->get(0),f);
}

void DL_multi_bar::begin_test(void) {
  DL_constraint::begin_test();
  for (int i=0;i<nr;i++)
    if (g_is_dyna[i]) ((DL_dyna *)g[i])->begintest();
}

void DL_multi_bar::end_test(void) {
  DL_constraint::end_test();
  for (int i=0;i<nr;i++)
    if (g_is_dyna[i]) ((DL_dyna *)g[i])->endtest();
}

boolean DL_multi_bar::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  static DL_largematrix dfdr(3,1);
  static DL_largematrix dcdf;
  static DL_largematrix subtemp;
  DL_vector ddiff;
  boolean nonzero=FALSE;
  dcdf.resize(cc->dim,3);
  subtemp.resize(cc->dim,dim);
//  sub->makezero();
  for (int i=0;i<nr;i++) {
    if (g_is_dyna[i]) {
      if (cc->dCdFq((DL_dyna*)g[i],&(p[i]),&dcdf)) {
        // calculate df/dr
        if (i==0) dfdr.setcolumn(0,&(d[0]));
	else
	if (i==nr-1) { dfdr.setcolumn(0,&(d[nr-2])); dfdr.neg(); }
	else {
          d[i].minus(&(d[i-1]),&ddiff);
	  dfdr.setcolumn(0,&ddiff);
	}
	// combine dc/df and df/dr into dc/dr:
        if (nonzero) {
	  dcdf.times(&dfdr,&subtemp);
	  sub->plusis(&subtemp);
	}
	else {
	  dcdf.times(&dfdr,sub);
	  nonzero=TRUE;
	}
      }
    }
  }
  return nonzero;
}

boolean DL_multi_bar::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  static DL_largematrix dcdp(1,3);
  static DL_largematrix dpdFq_(3,3);
  static DL_largematrix dcdf_term(1,3);
  DL_matrix dpdFq;
  DL_matrix ddpdFq;
  DL_vector ddiff;
  boolean nonzero=FALSE;
  for (int i=0; i<nr; i++) {
    if (g[i]==dc) {
      // calculate dp/df:
      dc->dpdfq(&(p[i]),pc,&dpdFq);
      if (veloterms) {
        dc->ddpdfq(&(p[i]),pc,&ddpdFq);
        ddpdFq.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdFq);
        dpdFq.plus(&ddpdFq,&dpdFq);
      }
      dpdFq_.assign(&dpdFq);
      // calculate dc/dp:
      if (i==0) { dcdp.setrow(0,&(d[0])); dcdp.neg(); }
      else
      if (i==nr-1) dcdp.setrow(0,&(d[nr-2]));
      else {
	d[i-1].minus(&(d[i]),&ddiff);
	dcdp.setrow(0,&ddiff);
      }
      // combine dc/dp and dp/df into dc/df
      if (nonzero) {
	dcdp.times(&dpdFq_,&dcdf_term);
	dcdfq->plusis(&dcdf_term);
      }
      else {
	dcdp.times(&dpdFq_,dcdfq);
	nonzero=TRUE;
      }
    }
  }
  return nonzero;
}

boolean DL_multi_bar::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  static DL_largematrix dcdp(1,3);
  static DL_largematrix dpdF_(3,3);
  static DL_largematrix dcdf_term(1,3);
  DL_matrix dpdF;
  DL_matrix ddpdF;
  DL_vector ddiff;
  boolean nonzero=FALSE;
  for (int i=0; i<nr; i++) {
    if (g[i]==dc) {
      // calculate dp/df:
      dc->dpdF(&(p[i]),&dpdF);
      if (veloterms) {
        dc->ddpdF(&(p[i]),&ddpdF);
        ddpdF.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdF);
        dpdF.plus(&ddpdF,&dpdF);
      }
      dpdF_.assign(&dpdF);
      // calculate dc/dp:
      if (i==0) { dcdp.setrow(0,&(d[0])); dcdp.neg(); }
      else
      if (i==nr-1) dcdp.setrow(0,&(d[nr-2]));
      else {
	d[i-1].minus(&(d[i]),&ddiff);
	dcdp.setrow(0,&ddiff);
      }
      // combine dc/dp and dp/df into dc/df
      if (nonzero) {
	dcdp.times(&dpdF_,&dcdf_term);
	dcdf->plusis(&dcdf_term);
      }
      else {
	dcdp.times(&dpdF_,dcdf);
	nonzero=TRUE;
      }
    }
  }
  return nonzero;
}

boolean DL_multi_bar::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  static DL_largematrix dcdp(1,3);
  static DL_largematrix dpdm_(3,3);
  static DL_largematrix dcdm_term(1,3);
  DL_matrix dpdm;
  DL_matrix ddpdm;
  DL_vector ddiff;
  boolean nonzero=FALSE;
  for (int i=0; i<nr; i++) {
    if (g[i]==dc) {
      // calculate dp/df:
      dc->dpdM(&(p[i]),&dpdm);
      if (veloterms) {
        dc->ddpdM(&(p[i]),&ddpdm);
        ddpdm.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdm);
        dpdm.plus(&ddpdm,&dpdm);
      }
      dpdm_.assign(&dpdm);
      // calculate dc/dp:
      if (i==0) { dcdp.setrow(0,&(d[0])); dcdp.neg(); }
      else
      if (i==nr-1) dcdp.setrow(0,&(d[nr-2]));
      else {
	d[i-1].minus(&(d[i]),&ddiff);
	dcdp.setrow(0,&ddiff);
      }
      // combine dc/dp and dp/dm into dc/dm
      if (nonzero) {
	 dcdp.times(&dpdm_,&dcdm_term);
	 dcdm->plusis(&dcdm_term);
      }
      else {
	dcdp.times(&dpdm_,dcdm);
	nonzero=TRUE;
      }
    }
  }
  return nonzero;
}

boolean DL_multi_bar::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  static DL_largematrix dcdp(1,3);
  static DL_largematrix dpdi_(3,3);
  static DL_largematrix dcdi_term(1,3);
  DL_matrix dpdi;
  DL_matrix ddpdi;
  DL_vector ddiff;
  boolean nonzero=FALSE;
  for (int i=0; i<nr; i++) {
    if (g[i]==dc) {
      // calculate dp/di:
      dc->dpdi(&(p[i]),pc,&dpdi);
      if (veloterms) {
        dc->ddpdi(&(p[i]),pc,&ddpdi);
        ddpdi.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdi);
        dpdi.plusis(&ddpdi);
      }
      dpdi_.assign(&dpdi);
      // calculate dc/dp:
      if (i==0) { dcdp.setrow(0,&(d[0])); dcdp.neg(); }
      else
      if (i==nr-1) dcdp.setrow(0,&(d[nr-2]));
      else {
	d[i-1].minus(&(d[i]),&ddiff);
	dcdp.setrow(0,&ddiff);
      }
      // combine dc/dp and dp/di into dc/di
      if (nonzero) {
	 dcdp.times(&dpdi_,&dcdi_term);
	 dcdi->plusis(&dcdi_term);
      }
      else {
	dcdp.times(&dpdi_,dcdi);
	nonzero=TRUE;
      }
    }
  }
  return nonzero;
}

DL_Scalar DL_multi_bar::actual_length(void){
  DL_vector pdiff;
  DL_point pdw,pgw;
  DL_Scalar len=0;
  if (g[0]) g[0]->to_world(&(p[0]),&pgw);
  else pgw.assign(&(p[0]));
  for(int i=1;i<nr;i++) {
    pdw.assign(&pgw);
    if (g[i]) g[i]->to_world(&(p[i]),&pgw);
    else pgw.assign(&(p[i]));
    pgw.minus(&pdw,&pdiff);
    len+=pdiff.norm();
  }
  return len;
}

void DL_multi_bar::test_restriction_changes(DL_largevector* lv) {
  if (maxforce>0) {
    if (fabs(F->get(0)+lv->get(0))>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: bar-constraint deactivated\n");
    }
  }
  if (rope) {
    if (F->get(0)+lv->get(0)<0) {
      // rope is trying to push: deactivate the bar and activate the
      // monitoring rope-controller:
      reset();
      rp.activate();
      deactivate();
    }
  }
}

void DL_multi_bar::new_frame(void) {
  DL_point pdw,pgw;
  int i;
  
  for (i=0; i<nr; i++) {
    if (g[i]) {
      g[i]->new_toworld(&(p[i]),&(pw[i]));
      g[i]->get_newvelocity(&(p[i]),&(dpw[i]));
    }
    else {
      p[i].minus(&(pw[i]),&(dpw[i]));
      dpw[i].timesis(1.0/DL_dsystem->get_integrator()->old_stepsize());
      pw[i].assign(&(p[i]));
    }
  }
  
  for (i=0;i<nr-1;i++) {
    pw[i+1].minus(&(pw[i]),&(d[i]));
    d[i].normalize();
  }

  DL_constraint::new_frame();
}

void DL_multi_bar::first_estimate(void) {
  static DL_largevector dF(dim);
  dF.resize(dim);
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
  
  if (check_restrictions()) apply_restrictions(F);
}

boolean DL_multi_bar::check_restrictions() {
  if (maxforce>0) {
    if (F->get(0)>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: bar-constraint deactivated\n");
      return FALSE;
    }
  }
  if (rope) {
    if (F->get(0)<0) {
      // rope is trying to push: deactivate the bar and activate the
      // monitoring rope-controller:
      reset();
      rp.activate();
      deactivate();
      return FALSE;
    }
  }
  return TRUE;
}

void DL_multi_bar::apply_restrictions(DL_largevector* lv) {
  DL_vector force;
  for (int i=0; i<nr-1; i++) {
    d[i].times(lv->get(0),&force);
    if (g_is_dyna[i]) ((DL_dyna*)g[i])->applyforce(&(p[i]),g[i],&force);
    if (g_is_dyna[i+1]) {
      force.neg(&force);
      ((DL_dyna*)g[i+1])->applyforce(&(p[i+1]),g[i+1],&force);
    }
  }
}

void DL_multi_bar::get_error(DL_largevector* lv) {
   DL_point pdw,pgw;
   DL_vector pdiff, dpdw, dpgw, vdiff;
   DL_Scalar err=-l;
   if (veloterms && (!testing)) err*=stiffness;

   if (g_is_dyna[0]) {
     g[0]->new_toworld(&(p[0]),&pgw);
     if (veloterms) g[0]->get_newvelocity(&(p[0]),&dpgw);
   }
   else {
     pgw.assign(&(pw[0]));
     if (veloterms) dpgw.assign(&(dpw[0]));
   }
   for (int i=1;i<nr;i++) {
     pdw.assign(&pgw);  
     if (g_is_dyna[i]) g[i]->new_toworld(&(p[i]),&pgw);
     else pgw.assign(&(pw[i]));
     pdw.minus(&pgw,&pdiff);
	 if (!testing) pdiff.timesis(stiffness);

     if (veloterms) {
       dpdw.assign(&dpgw);
       if (g_is_dyna[i]) g[i]->get_newvelocity(&(p[i]),&dpgw);
       else dpgw.assign(&(dpw[i]));
       dpdw.minus(&dpgw,&vdiff);
       vdiff.timesis(DL_dsystem->get_integrator()->halfstepsize());
       pdiff.plusis(&vdiff);
     }
     err+=pdiff.norm();
   }
   lv->set(0,err);
}

void DL_multi_bar::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  if (active) { nrf=nr-1; tf=2*nrf; }
  else { nrf=tf=0; }
}

void DL_multi_bar::get_force_info(int i, DL_actuator_type& at,
			         DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (active) {
    if ((0<=i) && (i<nr-1)) {
      if (g_is_dyna[i]) {
	at=force;
	_g=(DL_dyna*)(g[i]);
	_p->assign(&p[i]);
	reactionforce(i,_v);
	return;
      }
    }
    i-=nr-1;
    if ((0<=i) && (i<nr-1)) {
      if (g_is_dyna[i+1]) {
	at=force;
	_g=(DL_dyna*)(g[i+1]);
	_p->assign(&p[i+1]);
	reactionforce(i,_v);
	_v->neg(_v);
	return;
      }
    }
  }
  at=none;
};
