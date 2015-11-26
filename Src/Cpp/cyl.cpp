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
// filename     : cyl.cpp
// description	: non-in;ine methods of class DL_cyl
// author	: Bart Barenbrug     August '96
//

#include "dyna_system.h"
#include "cyl.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_cyl::DL_cyl():DL_constraint(),
  dcdp0(2,3),dcdp1(2,3),df0dR(3,2),df1dR(3,2) {
  dim=4;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  d=NULL;
  g=NULL;
  maxforce=0.0;
}

DL_cyl::~DL_cyl() {
}

void DL_cyl::init(DL_dyna* _d, DL_point* _pd0, DL_point *_pd1,
                          DL_geo* _g, DL_point* _pg0, DL_point *_pg1) {
  if (_d==_g) {
    DL_dsystem->get_companion()->Msg("Error: cylinder_constraint::init: a cylinder constraint needs points from _different_ objects\n cylinder constraint not initialised\n");
    return;
  }
  DL_constraint::init();
  d=_d;
  pd0.assign(_pd0);
  pd1.assign(_pd1);
  g=_g;
  pg0.assign(_pg0);
  pg1.assign(_pg1);
  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;
  if (g) {
    g->new_toworld(&pg0,&pg0w);
    g->new_toworld(&pg1,&pg1w);
    g->get_newvelocity(&pg0,&dpg0w);
    g->get_newvelocity(&pg1,&dpg1w);
  }
  else {
    pg0w.assign(&pg0);
    pg1w.assign(&pg1);
    dpg0w.init(0,0,0);
    dpg1w.init(0,0,0);
  }
  F->makezero(); oldF->makezero();

  // setup the local coordinate system:

  pg1w.minus(&pg0w,&l); l.normalize();

  // x and y are two vectors each with length 1 and both perpendicular to
  // each other and l
  if ((l.x!=0) || (l.z!=0)) x.init(l.z,0,-l.x);
  else x.init(0,0,-l.y);  x.normalize();
  l.crossprod(&x,&y);     y.normalize();
  
  dcdp0.setrow(0,&x);  df0dR.setcolumn(0,&x);
  dcdp0.setrow(1,&y);  df0dR.setcolumn(1,&y);
  dcdp1.setrow(0,&x);  df1dR.setcolumn(0,&x);
  dcdp1.setrow(1,&y);  df1dR.setcolumn(1,&y);

  // check if the constraint is initially valid
  DL_largevector lv(4);
  get_error(&lv);
  if (lv.norm()!=0.0) {
    DL_dsystem->get_companion()->Msg("Warning: initially invalid cylinder constraint.\n Error: (%f,%f,%f,%f)\n",
				  lv.get(0), lv.get(1), lv.get(2), lv.get(3) );
  }
}

void DL_cyl::getforce0(DL_largevector *lv, DL_vector *f) {
  DL_vector vtmp;
  x.times(lv->get(0),f);
  y.times(lv->get(1),&vtmp);
  f->plusis(&vtmp);
}

void DL_cyl::getforce1(DL_largevector *lv, DL_vector *f) {
  DL_vector vtmp;
  x.times(lv->get(2),f);
  y.times(lv->get(3),&vtmp);
  f->plusis(&vtmp);
}

void DL_cyl::reactionforce0(DL_vector *f) {
  getforce0(F,f);
}

void DL_cyl::reactionforce1(DL_vector *f) {
  getforce1(F,f);
}

DL_Scalar DL_cyl::reactionforce() {
  return F->norm();
}

void DL_cyl::begin_test(void) {
  DL_constraint::begin_test();
  d->begintest();
  if (g_is_dyna) ((DL_dyna *)g)->begintest();
}

void DL_cyl::end_test(void) {
  DL_constraint::end_test();
  d->endtest();
  if (g_is_dyna) ((DL_dyna*)g)->endtest();
}

boolean DL_cyl::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  static DL_largematrix dcdf0;
  static DL_largematrix dcdf1;
  static DL_largematrix subtemp;
  static DL_largematrix tmp;
  boolean nonzero;
  dcdf0.resize(cc->dim,3);
  if (nonzero=cc->dCdFq(d,&pd0,&dcdf0)) {
    dcdf1.resize(cc->dim,3);
    tmp.resize(cc->dim,2);
    dcdf0.times(&df0dR,&tmp);
    sub->setsubmatrix(0,0,&tmp);
    cc->dCdFq(d,&pd1,&dcdf1);
    dcdf1.times(&df1dR,&tmp);
    sub->setsubmatrix(0,2,&tmp);
  }
  if (g_is_dyna) {
   if (nonzero) {
    if (cc->dCdFq((DL_dyna*)g,&pg0,&dcdf0)) {
      dcdf1.resize(cc->dim,3);
      subtemp.resize(cc->dim,dim);
      tmp.resize(cc->dim,2);
      dcdf0.times(&df0dR,&tmp);
      subtemp.setsubmatrix(0,0,&tmp);
      cc->dCdFq((DL_dyna*)g,&pg1,&dcdf1);
      dcdf1.times(&df1dR,&tmp);
      subtemp.setsubmatrix(0,2,&tmp);
      sub->minusis(&subtemp);
    }
   }
   else {
    if (nonzero=cc->dCdFq((DL_dyna*)g,&pg0,&dcdf0)) {
      dcdf1.resize(cc->dim,3);
      tmp.resize(cc->dim,2);
      dcdf0.times(&df0dR,&tmp);
      sub->setsubmatrix(0,0,&tmp);
      cc->dCdFq((DL_dyna*)g,&pg1,&dcdf1);
      dcdf1.times(&df1dR,&tmp);
      sub->setsubmatrix(0,2,&tmp);
      sub->neg();
    }
   }
  }
  return nonzero;
}

boolean DL_cyl::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  static DL_largematrix dpdf(3,3);
  static DL_largematrix dcIdf(2,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;

  if (dc==d) { // there is an effect on C through d
    dc->dpdfq(&pd0,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&pd0,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp0.times(&dpdf,&dcIdf);
    dcdfq->setsubmatrix(0,0,&dcIdf);
    
    dc->dpdfq(&pd1,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&pd1,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp1.times(&dpdf,&dcIdf);
    dcdfq->setsubmatrix(2,0,&dcIdf);
    
    return TRUE;
  }
  if (dc==g) { // there is an effect on C through g
    dc->dpdfq(&pg0,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&pg0,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp0.times(&dpdf,&dcIdf);
    dcdfq->setsubmatrix(0,0,&dcIdf);
    
    dc->dpdfq(&pg1,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&pg1,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp1.times(&dpdf,&dcIdf);
    dcdfq->setsubmatrix(2,0,&dcIdf);
    
    dcdfq->neg();    
    return TRUE;
  }
  return FALSE;
}

boolean DL_cyl::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  static DL_largematrix dpdf(3,3);
  static DL_largematrix dcIdf(2,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dc->dpdF(&pd0,&dpdX);
    if (veloterms) {
      dc->ddpdF(&pd0,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp0.times(&dpdf,&dcIdf);
    dcdf->setsubmatrix(0,0,&dcIdf);
    
    dc->dpdF(&pd1,&dpdX);
    if (veloterms) {
      dc->ddpdF(&pd1,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp1.times(&dpdf,&dcIdf);
    dcdf->setsubmatrix(2,0,&dcIdf);
    
    return TRUE;
  }
  if (dc==g) { // there is an effect on C through g
    dc->dpdF(&pg0,&dpdX);
    if (veloterms) {
      dc->ddpdF(&pg0,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp0.times(&dpdf,&dcIdf);
    dcdf->setsubmatrix(0,0,&dcIdf);
    
    dc->dpdF(&pg1,&dpdX);
    if (veloterms) {
      dc->ddpdF(&pg1,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp1.times(&dpdf,&dcIdf);
    dcdf->setsubmatrix(2,0,&dcIdf);
    
    dcdf->neg();    
    return TRUE;
  }
  return FALSE;
}

boolean DL_cyl::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dpdm==0))
// dcdm has dimensions dim x 3
  static DL_largematrix dpdm(3,3);
  static DL_largematrix dcIdm(2,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dc->dpdM(&pd0,&dpdX);
    if (veloterms) {
      dc->ddpdM(&pd0,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdp0.times(&dpdm,&dcIdm);
    dcdm->setsubmatrix(0,0,&dcIdm);
    
    dc->dpdM(&pd1,&dpdX);
    if (veloterms) {
      dc->ddpdM(&pd1,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdp1.times(&dpdm,&dcIdm);
    dcdm->setsubmatrix(2,0,&dcIdm);
    
    return TRUE;
  }
  if (dc==g) { // there is an effect on C through g
    dc->dpdM(&pg0,&dpdX);
    if (veloterms) {
      dc->ddpdM(&pg0,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdp0.times(&dpdm,&dcIdm);
    dcdm->setsubmatrix(0,0,&dcIdm);
    
    dc->dpdM(&pg1,&dpdX);
    if (veloterms) {
      dc->ddpdM(&pg1,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdp1.times(&dpdm,&dcIdm);
    dcdm->setsubmatrix(2,0,&dcIdm);
    
    dcdm->neg();    
    return TRUE;
  }
  return FALSE;
}

boolean DL_cyl::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  static DL_largematrix dpdi(3,3);
  static DL_largematrix dcIdi(2,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dc->dpdi(&pd0,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&pd0,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdp0.times(&dpdi,&dcIdi);
    dcdi->setsubmatrix(0,0,&dcIdi);
    
    dc->dpdi(&pd1,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&pd1,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdp1.times(&dpdi,&dcIdi);
    dcdi->setsubmatrix(2,0,&dcIdi);
    
    return TRUE;
  }
  if (dc==g) { // there is an effect on C through g
    dc->dpdi(&pg0,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&pg0,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdp0.times(&dpdi,&dcIdi);
    dcdi->setsubmatrix(0,0,&dcIdi);
    
    dc->dpdi(&pg1,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&pg1,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdp1.times(&dpdi,&dcIdi);
    dcdi->setsubmatrix(2,0,&dcIdi);
    
    dcdi->neg();    
    return TRUE;
  }
  return FALSE;
}

void DL_cyl::get_error(DL_largevector* lv) {
   DL_point pdw;
   DL_vector dpdw, pdiff, vdiff;
   
   d->new_toworld(&pd0,&pdw);
   if (g_is_dyna) g->new_toworld(&pg0,&pg0w);
   pdw.minus(&pg0w,&pdiff);
   if (!testing) pdiff.timesis(stiffness);
   if (veloterms) {
     d->get_newvelocity(&pd0,&dpdw);
     if (g_is_dyna) g->get_newvelocity(&pg0,&dpg0w);
     dpdw.minus(&dpg0w,&vdiff);
     vdiff.timesis(DL_dsystem->get_integrator()->halfstepsize());
     pdiff.plusis(&vdiff);
   }
   lv->set(0,pdiff.inprod(&x));
   lv->set(1,pdiff.inprod(&y));
   
   d->new_toworld(&pd1,&pdw);
   if (g_is_dyna) g->new_toworld(&pg1,&pg1w);
   pdw.minus(&pg1w,&pdiff);
   if (!testing) pdiff.timesis(stiffness);
   if (veloterms) {
     d->get_newvelocity(&pd1,&dpdw);
     if (g_is_dyna) g->get_newvelocity(&pg1,&dpg1w);
     dpdw.minus(&dpg1w,&vdiff);
     vdiff.timesis(DL_dsystem->get_integrator()->halfstepsize());
     pdiff.plusis(&vdiff);
   }
   lv->set(2,pdiff.inprod(&x));
   lv->set(3,pdiff.inprod(&y));
}

void DL_cyl::test_restriction_changes(DL_largevector* lv) {
  if (maxforce>0) {
    DL_largevector Fnew(dim);
    F->plus(lv,&Fnew);
    if (Fnew.norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: cylinder-constraint deactivated\n");
    }
  }
}

boolean DL_cyl::check_restrictions() {
  if (maxforce>0) {
    if (F->norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: cylinder-constraint deactivated\n");
      return FALSE;
    }
  }
  return TRUE;
}

void DL_cyl::apply_restrictions(DL_largevector* lv) {
  DL_vector force0,force1;
  getforce0(lv,&force0);
  getforce1(lv,&force1);
  d->applyforce(&pd0,d,&force0);
  d->applyforce(&pd1,d,&force1);
  if (g_is_dyna) {
    force0.neg(&force0);
    force1.neg(&force1);
    ((DL_dyna*)g)->applyforce(&pg0,g,&force0);
    ((DL_dyna*)g)->applyforce(&pg1,g,&force1);
  }
}

void DL_cyl::first_estimate(void) {
  if (!g_is_dyna) {
    if (g) {
      g->new_toworld(&pg0,&pg0w);
      g->new_toworld(&pg1,&pg1w);
      g->get_newvelocity(&pg0,&dpg0w);
      g->get_newvelocity(&pg1,&dpg1w);
    }
    else {
      DL_Scalar hinv=1.0/DL_dsystem->get_integrator()->old_stepsize();
      pg0.minus(&pg0w,&dpg0w);
      dpg0w.timesis(hinv);
      pg0w.assign(&pg0);
      pg1.minus(&pg1w,&dpg1w);
      dpg1w.timesis(hinv);
      pg1w.assign(&pg1);
    }
  }

  { DL_vector oldl; oldl.assign(&l);
   
    pg1w.minus(&pg0w,&l);
    l.normalize();
    rotatebase(&oldl,&l,&x,&y);
  }
  
  dcdp0.setrow(0,&x);  df0dR.setcolumn(0,&x);
  dcdp0.setrow(1,&y);  df0dR.setcolumn(1,&y);
  dcdp1.setrow(0,&x);  df1dR.setcolumn(0,&x);
  dcdp1.setrow(1,&y);  df1dR.setcolumn(1,&y);

  DL_constraint::first_estimate();
}

void DL_cyl::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=2;
  if (g_is_dyna) tf=4; else tf=2;
}

void DL_cyl::get_force_info(int i, DL_actuator_type& at,
			      DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    at=force;
    _g=d;
    _p->assign(&pd0);
    reactionforce0(_v);
    return;
  }
  if (i==1) {
    at=force;
    _g=d;
    _p->assign(&pd1);
    reactionforce1(_v);
    return;
  }
  if (g_is_dyna) {
    if (i==2) {
      at=force;
      _g=(DL_dyna*)g;
      _p->assign(&pg0);
      reactionforce0(_v);
      _v->neg(_v);
      return;
    }
    if (i==3) {
      at=force;
      _g=(DL_dyna*)g;
      _p->assign(&pg1);
      reactionforce1(_v);
      _v->neg(_v);
      return;
    }
  }
  at=none;
};
