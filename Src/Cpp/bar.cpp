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
// filename     : ptp.cpp
// description	: non-inline methods of class DL_ptp
// author	: Bart Barenbrug     March '97
//

#include "dyna_system.h"
#include "bar.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_bar::DL_bar():DL_constraint(),
  dcdp(1,3), dfdr(3,1) {
  dim=1;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  maxforce=0.0;
  d=NULL;
  g=NULL;
  l=lsqr=1.0;
  rope=FALSE;
}

DL_bar::~DL_bar() {
}

void DL_bar::init(DL_dyna* _d, DL_point* _pd,
                          DL_geo* _g, DL_point* _pg,
			  DL_Scalar _l) {
  if (_d==_g) {
    DL_dsystem->get_companion()->Msg("Error: bar::init: a bar-constraint needs points from _different_ objects\n bar constraint not initialised\n");
    return;
  }
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Error: bar::init: for bars of length zero a ptp constraint should be used:\n bar-constraint not initialised\n");
    return;
  }
  DL_constraint::init();
  d=_d;
  pd.assign(_pd);
  g=_g;
  pg.assign(_pg);
  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;
  if (!g_is_dyna) {
    if (g) {
      g->new_toworld(&pg,&pgw);
      g->get_newvelocity(&pg,&dpgw);
    }
    else {
      pgw.assign(&pg);
      dpgw.init(0,0,0);
    }
  }
  l=_l;
  lsqr=l*l;
  F->makezero(); oldF->makezero();

  rp.init(this,&pd,&pg,lsqr);

  // check if the constraint is initially valid
  DL_largevector lv(1);
  get_error(&lv);
  if (rope) {
    if (lv.get(0)<0) {
      rp.activate();
      deactivate();
      return;
    }
  }
  if (lv.norm()>DL_constraints->max_error) {
    DL_dsystem->get_companion()->Msg("Warning: initially invalid bar-constraint.\n Error: %f\n", lv.get(0) );
  }

  DL_point pdw;
  DL_vector pdiff;
  d->to_world(&pd,&pdw);
  pdw.minus(&pgw,&pdiff);
  pdiff.timesis(2);
  dcdp.setrow(0,&pdiff);
  pdiff.normalize();
  dfdr.setcolumn(0,&pdiff);
}

void DL_bar::reactionforce(DL_vector *f) {
  DL_point pdw;
  DL_vector pdiff;
  
  d->to_world(&pd,&pdw);
  pdw.minus(&pgw,&pdiff);
  pdiff.normalize();
  pdiff.times(F->get(0),f);
}

void DL_bar::set_length(DL_Scalar _l){
  if (_l>0) {
    l=_l;
    lsqr=l*l;
    rp.set_lsqr(lsqr);
  }
  else {
    DL_dsystem->get_companion()->Msg("Warning: bar::set_length: length should be greater than zero! Old length kept\n");
  }
}

void DL_bar::begin_test(void) {
  DL_constraint::begin_test();
  d->begintest();
  if (g_is_dyna) ((DL_dyna *)g)->begintest();
}

void DL_bar::end_test(void) {
  DL_constraint::end_test();
  d->endtest();
  if (g_is_dyna) ((DL_dyna*)g)->endtest();
}

boolean DL_bar::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  static DL_largematrix dcdf;
  static DL_largematrix subtemp;
  boolean nonzero;
  dcdf.resize(cc->dim,3);
  if (nonzero=cc->dCdFq(d,&pd,&dcdf)) dcdf.times(&dfdr,sub);
  if (g_is_dyna) {
   if (nonzero) {
    if (cc->dCdFq((DL_dyna*)g,&pg,&dcdf)) {
      subtemp.resize(cc->dim,dim);
      dcdf.times(&dfdr,&subtemp);
      sub->minusis(&subtemp);
    }
   }
   else {
    if (nonzero=cc->dCdFq((DL_dyna*)g,&pg,&dcdf)) {
      dcdf.times(&dfdr,sub);
      sub->neg();
    }
   }
  }
  return nonzero;
}

boolean DL_bar::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  static DL_largematrix dpdFq_(3,3);
  DL_matrix dpdFq;
  DL_matrix ddpdFq;
  if (d==dc) {
    dc->dpdfq(&pd,pc,&dpdFq);
    if (veloterms) {
      dc->ddpdfq(&pd,pc,&ddpdFq);
      ddpdFq.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdFq);
      dpdFq.plus(&ddpdFq,&dpdFq);
    }
    dpdFq_.assign(&dpdFq);
    dcdp.times(&dpdFq_,dcdfq);
    return TRUE;
  }
  if (g==dc) {
    dc->dpdfq(&pg,pc,&dpdFq);
    if (veloterms) {
      dc->ddpdfq(&pg,pc,&ddpdFq);
      ddpdFq.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdFq);
      dpdFq.plus(&ddpdFq,&dpdFq);
    }
    dpdFq_.assign(&dpdFq);
    dcdp.times(&dpdFq_,dcdfq);
    dcdfq->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_bar::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  static DL_largematrix dpdf_(3,3);
  DL_matrix dpdf;
  DL_matrix ddpdf;
  if (d==dc) {
    dc->dpdF(&pd,&dpdf);
    if (veloterms) {
      dc->ddpdF(&pd,&ddpdf);
      ddpdf.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdf);
      dpdf.plus(&ddpdf,&dpdf);
    }
    dpdf_.assign(&dpdf);
    dcdp.times(&dpdf_,dcdf);
    return TRUE;
  }
  if (g==dc) {
    dc->dpdF(&pg,&dpdf);
    if (veloterms) {
      dc->ddpdF(&pg,&ddpdf);
      ddpdf.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdf);
      dpdf.plus(&ddpdf,&dpdf);
    }
    dpdf_.assign(&dpdf);
    dcdp.times(&dpdf_,dcdf);
    dcdf->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_bar::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  static DL_largematrix dpdm_(3,3);
  DL_matrix dpdm;
  DL_matrix ddpdm;
  if (d==dc) {
    dc->dpdM(&pd,&dpdm);
    if (veloterms) {
      dc->ddpdM(&pd,&ddpdm);
      ddpdm.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdm);
      dpdm.plus(&ddpdm,&dpdm);
    }
    dpdm_.assign(&dpdm);
    dcdp.times(&dpdm_,dcdm);
    return TRUE;
  }
  if (g==dc) {
    dc->dpdM(&pg,&dpdm);
    if (veloterms) {
      dc->ddpdM(&pg,&ddpdm);
      ddpdm.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdm);
      dpdm.plus(&ddpdm,&dpdm);
    }
    dpdm_.assign(&dpdm);
    dcdp.times(&dpdm_,dcdm);
    dcdm->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_bar::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  static DL_largematrix dpdi_(3,3);
  DL_matrix dpdi;
  DL_matrix ddpdi;
  if (d==dc) {
    dc->dpdi(&pd,pc,&dpdi);
    if (veloterms) {
      dc->ddpdi(&pd,pc,&ddpdi);
      ddpdi.timesis(DL_dsystem->get_integrator()->halfstepsize());
      dpdi.plus(&ddpdi,&dpdi);
    }
    dpdi_.assign(&dpdi);
    dcdp.times(&dpdi_,dcdi);
    return TRUE;
  }
  if (g==dc) {
    dc->dpdi(&pg,pc,&dpdi);
    if (veloterms) {
      dc->ddpdi(&pg,pc,&ddpdi);
      ddpdi.timesis(DL_dsystem->get_integrator()->halfstepsize());
      dpdi.plus(&ddpdi,&dpdi);
    }
    dpdi_.assign(&dpdi);
    dcdp.times(&dpdi_,dcdi);
    dcdi->neg();
    return TRUE;
  }
  return FALSE;
}

void DL_bar::new_frame(void) {
  if (!g_is_dyna) {
    if (g) {
      g->new_toworld(&pg,&pgw);
      g->get_newvelocity(&pg,&dpgw);
    }
    else {
      pg.minus(&pgw,&dpgw);
      dpgw.timesis(1.0/DL_dsystem->get_integrator()->old_stepsize());
      pgw.assign(&pg);
    }
  }
  
  DL_point pdw;
  DL_vector pdiff;
  d->to_world(&pd,&pdw);
  pdw.minus(&pgw,&pdiff);
  pdiff.timesis(2);
  dcdp.setrow(0,&pdiff);
  pdiff.normalize();
  dfdr.setcolumn(0,&pdiff);

  DL_constraint::new_frame();
}

void DL_bar::first_estimate(void) {
  if (oldF->norm()==0.0) {
    // constant extrapolation:
    oldF->assign(F);
  }
  else {
    // linear extrapolation:
    DL_largevector dF(dim);
    F->minus(oldF,&dF);
    oldF->assign(F);
    dF.timesis(stiffness);
    F->plusis(&dF);
  }
  
  if (check_restrictions()) apply_restrictions(F);
}

void DL_bar::test_restriction_changes(DL_largevector* lv) {
  if (maxforce>0) {
    if (fabs(F->get(0)+lv->get(0))>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: bar constraint deactivated\n");
    }
  }
  if (rope) {
    if ((F->get(0)+lv->get(0))>0) {
      // rope is trying to push: deactivate the bar and activate the
      // monitoring rope-controller:
      reset();
      rp.activate();
      deactivate();
    }
  }
}

boolean DL_bar::check_restrictions() {
  if (maxforce>0) {
    if (F->get(0)>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: bar-constraint deactivated\n");
      return FALSE;
    }
  }
  if (rope) {
    if (F->get(0)>0) {
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

void DL_bar::apply_restrictions(DL_largevector* lv) {
  DL_vector force;
  force.init(lv->get(0)*dfdr.get(0,0),
             lv->get(0)*dfdr.get(0,1),
	     lv->get(0)*dfdr.get(0,2));
  d->applyforce(&pd,d,&force);
  if (g_is_dyna) {
    force.neg(&force);
    ((DL_dyna*)g)->applyforce(&pg,g,&force);
  }
}

void DL_bar::get_error(DL_largevector* lv) {
   DL_point pdw;
   DL_vector pdiff;

   d->new_toworld(&pd,&pdw);
   if (g_is_dyna) g->new_toworld(&pg,&pgw);
   pdw.minus(&pgw,&pdiff);
   if (!testing) pdiff.timesis(stiffness);
   if (veloterms) {
     DL_vector dpdw, vdiff;
     d->get_newvelocity(&pd,&dpdw);
     if (g_is_dyna) g->get_newvelocity(&pg,&dpgw);
     dpdw.minus(&dpgw,&vdiff);
     vdiff.timesis(DL_dsystem->get_integrator()->halfstepsize());
     pdiff.plusis(&vdiff);
   }
   lv->init(pdiff.inprod(&pdiff)-lsqr);
}

void DL_bar::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=1;
  if (g_is_dyna) tf=2; else tf=1;
}

void DL_bar::get_force_info(int i, DL_actuator_type& at,
			    DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    at=force;
    _g=d;
    _p->assign(&pd);
    reactionforce(_v);
    return;
  }
  if ((i==1) && (g_is_dyna)) {
    at=force;
    _g=(DL_dyna*)g;
    _p->assign(&pg);
    reactionforce(_v);
    _v->neg(_v);
    return;
  }
  at=none;
};
