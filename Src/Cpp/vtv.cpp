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
// filename     : vtv.cpp
// description	: non-inline methods of class vtv
// author	: Bart Barenbrug     April '96
//

#include "dyna_system.h"
#include "vtv.h"

// ********************** //
// non-inline member fuctions //
// ********************** //

DL_vtv::DL_vtv():DL_constraint() {
  dim=3;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  maxforce=0.0;
  d=NULL;
  g=NULL;
}

DL_vtv::~DL_vtv() {
}

void DL_vtv::init(DL_dyna* _d, DL_point* _pd, DL_geo* _g, DL_point* _pg) {
  if (_d==_g) {
    DL_dsystem->get_companion()->Msg("Error: vtv::init: a vtv-constraint needs points from _different_ objects\n vtv-constraint not initialised\n");
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
  F->makezero(); oldF->makezero();

  // check if the constraint is initially valid
  DL_largevector lv(3);
  get_error(&lv);
  if (lv.norm()>DL_constraints->max_error) {
    DL_dsystem->get_companion()->Msg("Warning: initially invalid vtv-constraint.\n  error: (%f,%f,%f)\n",
				  lv.get(0), lv.get(1), lv.get(2) );
  }
}

void DL_vtv::reactionforce(DL_vector *f) {
  f->init(F->get(0), F->get(1), F->get(2));
}

void DL_vtv::begin_test(void) {
  DL_constraint::begin_test();
  d->begintest();
  if (g_is_dyna) ((DL_dyna *)g)->begintest();
}

void DL_vtv::end_test(void) {
  DL_constraint::end_test();
  d->endtest();
  if (g_is_dyna) ((DL_dyna*)g)->endtest();
}

boolean DL_vtv::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  boolean nonzero=cc->dCdFq(d,&pd,sub);
  if (g_is_dyna) {
    if (nonzero) {
      DL_largematrix subtemp(cc->dim,3);
      if (cc->dCdFq((DL_dyna*)g,&pg,&subtemp)) sub->minusis(&subtemp);
    }
    else if (nonzero=cc->dCdFq((DL_dyna*)g,&pg,sub)) sub->neg();
  }
  return nonzero;
}

boolean DL_vtv::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  DL_matrix ddpdFq;
  if (d==dc) {
    dc->ddpdfq(&pd,pc,&ddpdFq);
    dcdfq->assign(&ddpdFq);
    return TRUE;
  }
  if (g==dc) {
    dc->ddpdfq(&pg,pc,&ddpdFq);
    dcdfq->assign(&ddpdFq);
    dcdfq->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_vtv::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  DL_matrix ddpdf;
  if (d==dc) {
    dc->ddpdF(&pd,&ddpdf);
    dcdf->assign(&ddpdf);
    return TRUE;
  }
  if (g==dc) {
    dc->ddpdF(&pg,&ddpdf);
    dcdf->assign(&ddpdf);
    dcdf->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_vtv::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  DL_matrix ddpdm;
  if (d==dc) {
    dc->ddpdM(&pd,&ddpdm);
    dcdm->assign(&ddpdm);
    return TRUE;
  }
  if (g==dc) {
    dc->ddpdM(&pg,&ddpdm);
    dcdm->assign(&ddpdm);
    dcdm->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_vtv::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  DL_matrix ddpdi;
  if (d==dc) {
    dc->ddpdi(&pd,pc,&ddpdi);
    dcdi->assign(&ddpdi);
    return TRUE;
  }
  if (g==dc) {
    dc->ddpdi(&pg,pc,&ddpdi);
    dcdi->assign(&ddpdi);
    dcdi->neg();
    return TRUE;
  }
  return FALSE;
}

void DL_vtv::new_frame(void) {
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
  DL_constraint::new_frame();
}

void DL_vtv::test_restriction_changes(DL_largevector* lv) {
  if (maxforce>0) {
    DL_largevector Fnew(dim);
    F->plus(lv,&Fnew);
    if (Fnew.norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: vtv-constraint deactivated\n");
    }
  }
}

boolean DL_vtv::check_restrictions() {
  if (maxforce>0) {
    if (F->norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: vtv-constraint deactivated\n");
      return FALSE;
    }
  }
  return TRUE;
}

void DL_vtv::apply_restrictions(DL_largevector* lv) {
  DL_vector force;
  force.init(lv->get(0), lv->get(1), lv->get(2));
  d->applyforce(&pd,d,&force);
  if (g_is_dyna) {
    force.neg(&force);
    ((DL_dyna*)g)->applyforce(&pg,g,&force);
  }
}

void DL_vtv::get_error(DL_largevector* lv) {
   DL_point pdw, pdiff;
   DL_vector dpdw, vdiff;

   d->get_newvelocity(&pd,&dpdw);
   if (g_is_dyna) g->get_newvelocity(&pg,&dpgw);
   dpdw.minus(&dpgw,&vdiff);
   if (!testing) vdiff.timesis(stiffness);
   lv->init(vdiff.x,vdiff.y,vdiff.z);
}

void DL_vtv::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=1;
  if (g_is_dyna) tf=2; else tf=1;
}

void DL_vtv::get_force_info(int i, DL_actuator_type& at,
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
