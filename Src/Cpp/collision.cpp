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
// filename     : collision.cpp
// description	: non-inline methods of class DL_collision
// author	: Bart Barenbrug     July '96
//

#include "collision.h"
#include "dyna_system.h"
#include "constraint_manager.h"
#include "force_drawer.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_collision::DL_collision(DL_geo *_g0, DL_point *_p0,
                           DL_geo *_g1, DL_point *_p1,
                           DL_vector *_n, int mode):
                   DL_constraint(),dcdX(1,3),dXdR(3,1) {
  dim=1;
  oldF->resize(dim); oldF->makezero();
  F->resize(dim); F->makezero();
  Fsave->resize(dim);
  g0=g1=NULL;
  init(_g0,_p0, _g1,_p1, _n, mode);
}

void DL_collision::init(DL_geo *_g0, DL_point *_p0,
                        DL_geo *_g1, DL_point *_p1,
			DL_vector *_n, int mode) {
  // mode: <1: autodetect for positional constraint (dim=1/2)
  // mode:  1: dim=1
  // mode: >1: dim=2
  DL_vector vn,pdifft,pdiffth;
  DL_point p0t,p1t;
  DL_Scalar el,pt,vt;
  g0=_g0;
  if (g0) g0_is_dyna=g0->is_dyna();
  else g0_is_dyna=FALSE;
  g1=_g1;
  if (g1) g1_is_dyna=g1->is_dyna();
  else g1_is_dyna=FALSE;
  if (!(g0_is_dyna || g1_is_dyna)) {
    DL_dsystem->get_companion()->Msg("error: a collision between two non-dyna's occurred\n collision not handled\n");
    return;
  }
  DL_constraint::init();
  DL_constraints->add_collision(this);
  p0.assign(_p0);
  p1.assign(_p1);
  n.assign(_n);
  n.normalize();
  if (g0) {
    el=g0->get_elasticity();
    g0->get_velocity(_p0,&v0);
    g0->to_world(_p0,&p0t);
    g0->new_toworld(_p0,&p0w);
  }
  else {
    el=0;
    v0.init(0,0,0);
    p0w.assign(_p0);
    p0t.assign(_p0);
  }
  if (g1) {
    el+=g1->get_elasticity();
    g1->get_velocity(_p1,&v1);
    g1->to_world(_p1,&p1t);
    g1->new_toworld(_p1,&p1w);
  }
  else {
    v1.init(0,0,0);
    p1w.assign(_p1);
    p1t.assign(_p1);
  }
  if (g0 && g1) el*=0.5;
  v1.minus(&v0,&vn);
  vt=vn.inprod(&n);
  v=el*vt;

  p1t.minus(&p0t,&pdifft);
  pt=pdifft.inprod(&n);
  p1w.minus(&p0w,&pdiffth);
  p=pdiffth.inprod(&n);
  pdiffth.minusis(&pdifft);
  if ((mode>=2) ||
      ((mode<=0) &&
       ((fabs(fabs(pt/(pt-p))-0.5)>0.25) || (el<0.5))
      )
     ){
    // add positional constraint.
    if (p*pt>0) {
      pt=p-2*pt;
      if (pt*p<0) p=0;
      else p=pt;
    }
    p*=el;
    dim=2;
  }
  else dim=1;
  // v0 and v1 were used for velocity at t until now: switch to t+h:
  if (g0 && !g0_is_dyna) g0->get_newvelocity(_p0,&v0);
  if (g1 && !g1_is_dyna) g1->get_newvelocity(_p1,&v1);
  dcdX.setrow(0,&n);
  dXdR.setcolumn(0,&n);
  oldF->resize(dim); oldF->makezero();
  F->resize(dim); F->makezero();
  Fsave->resize(dim);
}

void DL_collision::begin_test(void) {
  DL_constraint::begin_test();
  if (g0_is_dyna) ((DL_dyna*)g0)->begintest();
  if (g1_is_dyna) ((DL_dyna*)g1)->begintest();
}

void DL_collision::end_test(void) {
  DL_constraint::end_test();
  if (g0_is_dyna) ((DL_dyna*)g0)->endtest();
  if (g1_is_dyna) ((DL_dyna*)g1)->endtest();
}

boolean DL_collision::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  boolean nonzero=FALSE;
  DL_largematrix dcdx(cc->dim,3); dcdx.makezero();
  if (g0_is_dyna) nonzero=cc->dCdFq((DL_dyna*)g0,&p0,&dcdx);
  if (g1_is_dyna) {
    DL_largematrix dcdx1(cc->dim,3);
    if (cc->dCdFq((DL_dyna*)g1,&p1,&dcdx1)) {
      dcdx.minusis(&dcdx1);
      nonzero=TRUE;
    }
  }
  if (nonzero) {
    DL_largematrix dcdrsub(cc->dim,1);
    dcdx.times(&dXdR,&dcdrsub);
    sub->setsubmatrix(0,0,&dcdrsub);
  }
  if (dim==1) return nonzero;
  
  if (g0_is_dyna) {
    if (cc->dCdI((DL_dyna*)g0,&p0,&dcdx)) nonzero=TRUE;
  }
  if (g1_is_dyna) {
    DL_largematrix dcdx1(cc->dim,3);
    if (cc->dCdI((DL_dyna*)g1,&p1,&dcdx1)) {
      dcdx.minusis(&dcdx1);
      nonzero=TRUE;
    }
  }
  if (nonzero) {
    DL_largematrix dcdrsub(cc->dim,1);
    dcdx.times(&dXdR,&dcdrsub);
    sub->setsubmatrix(0,1,&dcdrsub);
  }
  return nonzero;
}

boolean DL_collision::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  if (dc==g0) {// which implies g0_is_dyna; there is an effect on C through g0
    DL_largematrix dvdf(3,3);
    DL_matrix ddpdf;
    dc->ddpdfq(&p0,pc,&ddpdf);
    dvdf.assign(&ddpdf);
    dvdf.neg();
    if (dim==1)
      dcdX.times(&dvdf,dcdfq);
    else {
      DL_largematrix dcdfqsub(1,3);
      dcdX.times(&dvdf,&dcdfqsub);
      dcdfq->setsubmatrix(0,0,&dcdfqsub);
      dc->dpdfq(&p0,pc,&ddpdf);
      dvdf.assign(&ddpdf);
      dvdf.neg();
      dcdX.times(&dvdf,&dcdfqsub);
      dcdfq->setsubmatrix(1,0,&dcdfqsub);
    }
    return TRUE;
  }
  if (dc==g1) {// which implies g1_is_dyna; there is an effect on C through g1
    DL_largematrix dvdf(3,3);
    DL_matrix ddpdf;
    dc->ddpdfq(&p1,pc,&ddpdf);
    dvdf.assign(&ddpdf);
    if (dim==1)
      dcdX.times(&dvdf,dcdfq);
    else {
      DL_largematrix dcdfqsub(1,3);
      dcdX.times(&dvdf,&dcdfqsub);
      dcdfq->setsubmatrix(0,0,&dcdfqsub);
      dc->dpdfq(&p1,pc,&ddpdf);
      dvdf.assign(&ddpdf);
      dcdX.times(&dvdf,&dcdfqsub);
      dcdfq->setsubmatrix(1,0,&dcdfqsub);
    }
    return TRUE;
  }
  return FALSE;
}

boolean DL_collision::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  if (dc==g0) {// which implies g0_is_dyna; there is an effect on C through g0
    DL_largematrix dvdf(3,3);
    DL_matrix ddpdf;
    dc->ddpdF(&p0,&ddpdf);
    dvdf.assign(&ddpdf);
    dvdf.neg();
    if (dim==1)
      dcdX.times(&dvdf,dcdf);
    else {
      DL_largematrix dcdfsub(1,3);
      dcdX.times(&dvdf,&dcdfsub);
      dcdf->setsubmatrix(0,0,&dcdfsub);
      dc->dpdF(&p0,&ddpdf);
      dvdf.assign(&ddpdf);
      dvdf.neg();
      dcdX.times(&dvdf,&dcdfsub);
      dcdf->setsubmatrix(1,0,&dcdfsub);
    }
    return TRUE;
  }
  if (dc==g1) {// which implies g1_is_dyna; there is an effect on C through g1
    DL_largematrix dvdf(3,3);
    DL_matrix ddpdf;
    dc->ddpdF(&p1,&ddpdf);
    dvdf.assign(&ddpdf);
    if (dim==1)
      dcdX.times(&dvdf,dcdf);
    else {
      DL_largematrix dcdfsub(1,3);
      dcdX.times(&dvdf,&dcdfsub);
      dcdf->setsubmatrix(0,0,&dcdfsub);
      dc->dpdF(&p1,&ddpdf);
      dvdf.assign(&ddpdf);
      dcdX.times(&dvdf,&dcdfsub);
      dcdf->setsubmatrix(1,0,&dcdfsub);
    }
    return TRUE;
  }
  return FALSE;
}

boolean DL_collision::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  if (dc==g0) {// which implies g0_is_dyna; there is an effect on C through g0
    DL_largematrix dvdm(3,3);
    DL_matrix ddpdm;
    dc->ddpdM(&p0,&ddpdm);
    dvdm.assign(&ddpdm);
    dvdm.neg();
    if (dim==1)
      dcdX.times(&dvdm,dcdm);
    else {
      DL_largematrix dcdmsub(1,3);
      dcdX.times(&dvdm,&dcdmsub);
      dcdm->setsubmatrix(0,0,&dcdmsub);
      dc->dpdM(&p0,&ddpdm);
      dvdm.assign(&ddpdm);
      dvdm.neg();
      dcdX.times(&dvdm,&dcdmsub);
      dcdm->setsubmatrix(1,0,&dcdmsub);
    }
    return TRUE;
  }
  if (dc==g1) {// which implies g1_is_dyna; there is an effect on C through g1
    DL_largematrix dvdm(3,3);
    DL_matrix ddpdm;
    dc->ddpdM(&p1,&ddpdm);
    dvdm.assign(&ddpdm);
    if (dim==1)
      dcdX.times(&dvdm,dcdm);
    else {
      DL_largematrix dcdmsub(1,3);
      dcdX.times(&dvdm,&dcdmsub);
      dcdm->setsubmatrix(0,0,&dcdmsub);
      dc->dpdM(&p1,&ddpdm);
      dvdm.assign(&ddpdm);
      dcdX.times(&dvdm,&dcdmsub);
      dcdm->setsubmatrix(1,0,&dcdmsub);
    }
    return TRUE;
  }
  return FALSE;
}

boolean DL_collision::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  if (dc==g0) {// which implies g0_is_dyna; there is an effect on C through g0
    DL_largematrix dvdI(3,3);
    DL_matrix ddpdi;
    dc->ddpdi(&p0,pc,&ddpdi);
    dvdI.assign(&ddpdi);
    dvdI.neg();
    if (dim==1)
      dcdX.times(&dvdI,dcdi);
    else {
      DL_largematrix dcdisub(1,3);
      dcdX.times(&dvdI,&dcdisub);
      dcdi->setsubmatrix(0,0,&dcdisub);
      dc->dpdi(&p0,pc,&ddpdi);
      dvdI.assign(&ddpdi);
      dvdI.neg();
      dcdX.times(&dvdI,&dcdisub);
      dcdi->setsubmatrix(1,0,&dcdisub);
    }
    return TRUE;
  }
  if (dc==g1) {// which implies g1_is_dyna; there is an effect on C through g1
    DL_largematrix dvdI(3,3);
    DL_matrix ddpdi;
    dc->ddpdi(&p1,pc,&ddpdi);
    dvdI.assign(&ddpdi);
    if (dim==1)
      dcdX.times(&dvdI,dcdi);
    else {
      DL_largematrix dcdisub(1,3);
      dcdX.times(&dvdI,&dcdisub);
      dcdi->setsubmatrix(0,0,&dcdisub);
      dc->dpdi(&p1,pc,&ddpdi);
      dvdI.assign(&ddpdi);
      dcdX.times(&dvdI,&dcdisub);
      dcdi->setsubmatrix(1,0,&dcdisub);
    }
    return TRUE;
  }
  return FALSE;
}

void DL_collision::apply_restrictions(DL_largevector* lv) {
  DL_vector force;
  
  n.times(lv->get(0),&force);
  if (g0_is_dyna) ((DL_dyna*)g0)->applyforce(&p0,g0,&force);
  if (g1_is_dyna) {
    force.neg(&force);
    ((DL_dyna*)g1)->applyforce(&p1,g1,&force);
  }
  
  if (dim==1) return;
  
  n.times(lv->get(1),&force); // use force to store the impulse now
  if (g0_is_dyna) ((DL_dyna*)g0)->applyimpulse(&p0,g0,&force);
  if (g1_is_dyna) {
    force.neg(&force);
    ((DL_dyna*)g1)->applyimpulse(&p1,g1,&force);
  }
}

void DL_collision::get_error(DL_largevector* lv) {
  DL_vector vdiff;
  if (g0_is_dyna) g0->get_newvelocity(&p0,&v0);
  if (g1_is_dyna) g1->get_newvelocity(&p1,&v1);
  v1.minus(&v0,&vdiff);
  lv->set(0,vdiff.inprod(&n)+v);

  if (dim==1) return;
  
  if (g0_is_dyna) g0->new_toworld(&p0,&p0w);    
  if (g1_is_dyna) g1->new_toworld(&p1,&p1w);
  p1w.minus(&p0w,&vdiff);
  lv->set(1,vdiff.inprod(&n)+p);
}

void DL_collision::post_processing(void) {
  // if the constraint manager hasn't satisfied the constraint with enough
  // accuracy, here is the place to test that and to perhaps create a
  // `follow-up' constraint for the next frame

  // show the collision force?
  if (DL_constraints->showing_constraint_forces()) {
    DL_vector force;
    n.times(F->get(0),&force);
    if (g0_is_dyna)
      DL_dsystem->get_companion()->draw_force((DL_dyna*)g0,&p0,&force);
    if (g1_is_dyna) {
      force.neg(&force);
      DL_dsystem->get_companion()->draw_reaction_force((DL_dyna*)g1,&p1,&force);
    }

    if (dim==2) {
      n.times(F->get(1),&force); // use "force" to store the impulse
      if (g0_is_dyna)
	DL_dsystem->get_companion()->draw_impulse((DL_dyna*)g0,&p0,&force);
      if (g1_is_dyna) {
	force.neg(&force);
	DL_dsystem->get_companion()->draw_reaction_impulse((DL_dyna*)g1,&p1,&force);
      }
    }
  }
  
  delete this;
}
