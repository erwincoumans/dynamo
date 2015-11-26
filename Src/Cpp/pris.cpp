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
// filename     : pris.cpp
// description	: non-inline methods of the DL_prism constraint class
// author	: Bart Barenbrug     August '96
//

#include "dyna_system.h"
#include "pris.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_pris::DL_pris():DL_constraint(),
  dcdp(2,3), dfdR(3,2) {
  dim=2;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  maxforce=0.0;
  d=NULL; g=NULL;
  g_is_dyna=FALSE;
  myorient=new DL_orientation();
}

DL_pris::~DL_pris() {
}

void DL_pris::init(DL_dyna *_d, DL_point *_pd, DL_vector *_ld, DL_vector *_rd,
                       DL_geo  *_g, DL_point* _pg, DL_vector *_lg, DL_vector* _rg) {
  if (_g==_d) {
    DL_dsystem->get_companion()->Msg("Error: prism_constraint::init: two different geometries are required!\n prism-constraint not initialised\n");
    return;
  }

  d=_d;
  g=_g;
  pd.assign(_pd);
  pg.assign(_pg);

  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;

  // set up the directions for the orientation constraint:
  DL_Scalar tf; DL_vector w0; tf=d->torquefactor();

  v0.assign(_ld); v0.normalize();
  v1.assign(_rd); v1.normalize();

  w0.assign(_lg); w0.normalize();
  w1.assign(_rg); w1.normalize();
  w0.crossprod(&w1,&w2);
  w2.normalize();
  
  myorient->init(_d,&v0,&v1,_g,&w1,&w2);  
  DL_constraint::init();
  F->makezero(); oldF->makezero();

  // set up the local coordinate system:
  d->to_world(&v0,&l);
  if (g) {
    g->new_toworld(&pg,&pgw);
    g->get_newvelocity(&pg,&dpgw);
    g->to_world(&w1,&x);
    g->to_world(&w2,&y);
  }
  else {
    pgw.assign(&pg);
    x.assign(&w1);
    y.assign(&w2);
    dpgw.init(0,0,0);
  }

  dcdp.setrow(0,&x);
  dcdp.setrow(1,&y);
  dfdR.setcolumn(0,&x);
  dfdR.setcolumn(1,&y);

  // check if the constraint is initially valid
  DL_largevector lv(2);
  get_error(&lv);
    
  if (lv.norm()>DL_constraints->max_error) {
    DL_dsystem->get_companion()->Msg("Warning: initially invalid prism-constraint\n Error: (%f,%f)\n",
				  lv.get(0), lv.get(1) );
  }
}

void DL_pris::getforce(DL_largevector* lv, DL_vector *force) {
  DL_vector vtmp;
  x.times(lv->get(0),force);
  y.times(lv->get(1),&vtmp);
  force->plusis(&vtmp);
}

void DL_pris::reactionforce(DL_vector *force) {
  getforce(F,force);
}

void DL_pris::begin_test(void) {
  DL_constraint::begin_test();
  d->begintest();
  if (g_is_dyna) ((DL_dyna *)g)->begintest();
}

void DL_pris::end_test(void) {
  DL_constraint::end_test();
  d->endtest();
  if (g_is_dyna) ((DL_dyna *)g)->endtest();
}

boolean DL_pris::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  static DL_largematrix dcdf;
  static DL_largematrix dcdfg;
  dcdf.resize(cc->dim,3);
  boolean nonzero=cc->dCdFq((DL_dyna*)d,&pd,&dcdf);
  if (g_is_dyna) {
    if (nonzero) {
      dcdfg.resize(cc->dim,3);
      if (cc->dCdFq((DL_dyna*)g,&pg,&dcdfg)) dcdf.minusis(&dcdfg);
    }
    else {
      if (nonzero=cc->dCdFq((DL_dyna*)g,&pg,&dcdf)) dcdf.neg();
    }
  }
  if (nonzero) dcdf.times(&dfdR,sub);
  return nonzero;
}

boolean DL_pris::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  static DL_largematrix dpdfq(3,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dc->dpdfq(&pd,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&pd,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdfq.assign(&dpdX);
    dcdp.times(&dpdfq,dcdfq);
    return TRUE;
  }
  if (dc==g) { // which implies g_is_dyna;
               // there is an effect on C through g
    dc->dpdfq(&pg,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&pg,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdfq.assign(&dpdX);
    dpdfq.neg();
    dcdp.times(&dpdfq,dcdfq);
    return TRUE;
  }
  return FALSE;
}

boolean DL_pris::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  static DL_largematrix dpdf(3,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dc->dpdF(&pd, &dpdX);
    if (veloterms) {
      dc->ddpdF(&pd, &ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp.times(&dpdf,dcdf);
    return TRUE;
  }
  if (dc==g) { // which implies g_is_dyna;
               // there is an effect on C through g
    dc->dpdF(&pg, &dpdX);
    if (veloterms) {
      dc->ddpdF(&pg, &ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dpdf.neg();
    dcdp.times(&dpdf,dcdf);
    return TRUE;
  }
  return FALSE;
}

boolean DL_pris::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  static DL_largematrix dpdm(3,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dc->dpdM(&pd, &dpdX);
    if (veloterms) {
      dc->ddpdM(&pd, &ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdp.times(&dpdm,dcdm);
    return TRUE;
  }
  if (dc==g) { // which implies g_is_dyna;
               // there is an effect on C through g
    dc->dpdM(&pg, &dpdX);
    if (veloterms) {
      dc->ddpdM(&pg, &ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dpdm.neg();
    dcdp.times(&dpdm,dcdm);
    return TRUE;
  }
  return FALSE;
}

boolean DL_pris::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  static DL_largematrix dpdi(3,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dc->dpdi(&pd,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&pd,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdp.times(&dpdi,dcdi);
    return TRUE;
  }
  if (dc==g) { // which implies g_is_dyna;
               // there is an effect on C through g
    dc->dpdi(&pg,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&pg,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dpdi.neg();
    dcdp.times(&dpdi,dcdi);
    return TRUE;
  }
  return FALSE;
}

void DL_pris::new_frame(void) {
   d->to_world(&v0,&l);
   if (g) {
     g->to_world(&w1,&x);
     g->to_world(&w2,&y);
   }

   dcdp.setrow(0,&x);
   dcdp.setrow(1,&y);
   dfdR.setcolumn(0,&x);
   dfdR.setcolumn(1,&y);

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

void DL_pris::test_restriction_changes(DL_largevector* lv) {
  if (!myorient->active) {
    deactivate();
    // possibly raise an event here
    DL_dsystem->get_companion()->Msg(" prism constraint deactivated\n");
  }
  if (maxforce>0) {
    DL_largevector Fnew(dim);
    F->plus(lv,&Fnew);
    if (Fnew.norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: prism-constraint deactivated\n");
    }
  }
}

boolean DL_pris::check_restrictions() {
  if (!myorient->active) {
    deactivate();
    // possibly raise an event here
    DL_dsystem->get_companion()->Msg(" prism constraint deactivated\n");
    return FALSE;
  }
  if (maxforce>0)
    if (F->norm()>maxforce) {
      myorient->deactivate();
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: prism-constraint deactivated\n");
      return FALSE;
    }
  return TRUE;
}

void DL_pris::apply_restrictions(DL_largevector* lv) {
  DL_vector force;
  getforce(lv,&force);
  d->applyforce(&pd,d,&force);
  if (g_is_dyna) {
    force.neg(&force);
    ((DL_dyna*)g)->applyforce(&pg,g,&force);
  }
}

void DL_pris::get_error(DL_largevector* lv) {

   DL_point pdw;
   DL_vector pdiff, vdiff;
   
   d->new_toworld(&pd,&pdw);
   if (g_is_dyna) g->new_toworld(&pg,&pgw);
   pdw.minus(&pgw,&pdiff);
   if (!testing) pdiff.timesis(stiffness);
   if (veloterms) {
     DL_vector dpdw;
     d->get_newvelocity(&pd,&dpdw);
     if (g_is_dyna) g->get_newvelocity(&pg,&dpgw);
     dpdw.minus(&dpgw,&vdiff);
     vdiff.timesis(DL_dsystem->get_integrator()->halfstepsize());
     pdiff.plusis(&vdiff);
   }
   lv->init(pdiff.inprod(&x), pdiff.inprod(&y));
}

void DL_pris::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=1;
  if (g_is_dyna) tf=2; else tf=1;
}

void DL_pris::get_force_info(int i, DL_actuator_type& at,
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
