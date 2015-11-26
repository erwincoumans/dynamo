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
// filename     : plc.cpp
// description	: non-inline methods of class DL_plc
// author	: Bart Barenbrug     August 1996
//

#include "dyna_system.h"
#include "plc.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_plc::DL_plc():DL_constraint(),
  dcdp(1,3), dcdv0(2,3), dcdw(1,3), dfdR(3,1), dmdR(3,2) {
  dim=3;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  d=NULL; g=NULL;
  maxforce=maxtorque=0.0;
}

DL_plc::~DL_plc() {
}

void DL_plc::init(DL_dyna* _d, DL_point* _pd, DL_vector *_ld,
                  DL_geo* _g,  DL_point* _pg, DL_vector *_lg) {
  if (_d==_g) {
    DL_dsystem->get_companion()->Msg("Error: plane_constraint::init: a plane constraint needs points from _different_ objects\n plane-constraint not initialised\n");
    return;
  }
  DL_constraint::init();
  d=_d;
  pd.assign(_pd);
  g=_g;
  pg.assign(_pg);

  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;

  DL_vector w0; DL_Scalar lngth,tf=d->torquefactor();

  v0.assign(_ld);    
  if ((lngth=v0.norm())!=0) v0.timesis(1.0/lngth);

  w0.assign(_lg);
  if ((lngth=w0.norm())!=0) w0.timesis(1.0/lngth);
  if ((w0.x!=0) || (w0.z!=0)) w1.init(w0.z,0,-w0.x);
  else w1.init(0,0,-w0.y);
  if ((lngth=w1.norm())!=0) w1.timesis(tf/lngth);
  w0.crossprod(&w1,&w2);
  if ((lngth=w2.norm())!=0) w2.timesis(tf/lngth);

  d->to_world(&v0,&l);
  if (g) {
    g->new_toworld(&pg,&pgw);
    g->get_newvelocity(&pg,&dpgw);
    g->new_toworld(&w1,&y);
    g->new_toworld(&w2,&x);
    g->get_newvelocity(&w1,&dw1w);
    g->get_newvelocity(&w2,&dw2w);
  }
  else {
    pgw.assign(&pg);
    dpgw.init(0,0,0);
    y.assign(&w1);
    x.assign(&w2);
    w1w.assign(&x);
    w2w.assign(&y);
    dw1w.init(0,0,0);
    dw2w.init(0,0,0);
  }

  F->makezero(); oldF->makezero();

  // check if the constraint is initially valid
  DL_largevector lv(3);
  get_error(&lv);
  if (lv.norm()>DL_constraints->max_error) {
    DL_dsystem->get_companion()->Msg("Warning: initially invalid plane-constraint.\n Error: (%f,%f,%f)\n",
				     lv.get(0), lv.get(1), lv.get(2) );
  }

  dfdR.setcolumn(0,&l);
  dmdR.setcolumn(0,&x);
  dmdR.setcolumn(1,&y);
  dcdp.setrow(0,&l);
  dcdv0.setrow(0,&y);
  dcdv0.setrow(1,&x);
  dcdw.setrow(0,&l);
}

void DL_plc::getforce(DL_largevector *lv, DL_vector *f) {
  l.times(lv->get(0),f);
}

void DL_plc::gettorque(DL_largevector *lv, DL_vector *t) {
  DL_vector vtmp;
  x.times(lv->get(1),t);
  y.times(lv->get(2),&vtmp);
  t->plusis(&vtmp);
}

void DL_plc::reactionforce(DL_vector *f) {
  getforce(F,f);
}

void DL_plc::reactiontorque(DL_vector *t) {
  gettorque(F,t);
}

void DL_plc::begin_test(void) {
  DL_constraint::begin_test();
  d->begintest();
  if (g_is_dyna) ((DL_dyna *)g)->begintest();
}

void DL_plc::end_test(void) {
  DL_constraint::end_test();
  d->endtest();
  if (g_is_dyna) ((DL_dyna*)g)->endtest();
}

boolean DL_plc::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  static DL_largematrix dcdf;
  static DL_largematrix subtemp;
  static DL_largematrix tmp;
  boolean nonzero;
  dcdf.resize(cc->dim,3);
  if (nonzero=cc->dCdFq(d,&pd,&dcdf)) {
    tmp.resize(cc->dim,1);
    dcdf.times(&dfdR,&tmp);
    sub->setsubmatrix(0,0,&tmp);
    tmp.resize(cc->dim,2);
    cc->dCdM(d,&dcdf);
    dcdf.times(&dmdR,&tmp);
    sub->setsubmatrix(0,1,&tmp);
  }
  if (g_is_dyna) {
   if (nonzero) {
    if (cc->dCdFq((DL_dyna*)g,&pg,&dcdf)) {
      subtemp.resize(cc->dim,dim);
      tmp.resize(cc->dim,1);
      dcdf.times(&dfdR,&tmp);
      subtemp.setsubmatrix(0,0,&tmp);
      tmp.resize(cc->dim,2);
      cc->dCdM((DL_dyna*)g,&dcdf);
      dcdf.times(&dmdR,&tmp);
      subtemp.setsubmatrix(0,1,&tmp);
      sub->minusis(&subtemp);
    }
   }
   else {
    if (nonzero=cc->dCdFq((DL_dyna*)g,&pg,&dcdf)) {
      tmp.resize(cc->dim,1);
      dcdf.times(&dfdR,&tmp);
      sub->setsubmatrix(0,0,&tmp);
      tmp.resize(cc->dim,2);
      cc->dCdM((DL_dyna*)g,&dcdf);
      dcdf.times(&dmdR,&tmp);
      sub->setsubmatrix(0,1,&tmp);
      sub->neg();
    }
   }
  }
  return nonzero;
}

boolean DL_plc::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  static DL_largematrix dpdf(3,3);
  static DL_largematrix dcIdf(1,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dcIdf.resize(1,3);
    dc->dpdfq(&pd,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&pd,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp.times(&dpdf,&dcIdf);
    dcdfq->setsubmatrix(0,0,&dcIdf);
    
    dc->dvdfq(&v0,pc,&dpdX);
    if (veloterms) {
      dc->ddvdfq(&v0,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcIdf.resize(2,3);
    dcdv0.times(&dpdf,&dcIdf);
    dcdfq->setsubmatrix(1,0,&dcIdf);
    
    return TRUE;
  }
  if (dc==g) { // there is an effect on C through g
    dcIdf.resize(1,3);
    dc->dpdfq(&pg,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&pg,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp.times(&dpdf,&dcIdf);
    dcIdf.neg();
    dcdfq->setsubmatrix(0,0,&dcIdf);
    
    dc->dvdfq(&w1,pc,&dpdX);
    if (veloterms) {
      dc->ddvdfq(&w1,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdw.times(&dpdf,&dcIdf);
    dcdfq->setsubmatrix(1,0,&dcIdf);
    
    dc->dvdfq(&w2,pc,&dpdX);
    if (veloterms) {
      dc->ddvdfq(&w2,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdw.times(&dpdf,&dcIdf);
    dcdfq->setsubmatrix(2,0,&dcIdf);
    
    return TRUE;
  }
  return FALSE;
}

boolean DL_plc::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  static DL_largematrix dpdf(3,3);
  static DL_largematrix dcIdf(1,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    // dcIdf.resize(1,3);
    dc->dpdF(&pd,&dpdX);
    if (veloterms) {
      dc->ddpdF(&pd,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp.times(&dpdf,&dcIdf);
    dcdf->setsubmatrix(0,0,&dcIdf);
    
    dcIdf.makezero(); // no effect of central force on orientation
    dcdf->setsubmatrix(1,0,&dcIdf);
    dcdf->setsubmatrix(2,0,&dcIdf);
    
    return TRUE;
  }
  if (dc==g) { // there is an effect on C through g
    // dcIdf.resize(1,3);
    dc->dpdF(&pg,&dpdX);
    if (veloterms) {
      dc->ddpdF(&pg,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp.times(&dpdf,&dcIdf);
    dcIdf.neg();
    dcdf->setsubmatrix(0,0,&dcIdf);
    
    dcIdf.makezero(); // no effect of central force on orientation
    dcdf->setsubmatrix(1,0,&dcIdf);
    dcdf->setsubmatrix(2,0,&dcIdf);
    
    return TRUE;
  }
  return FALSE;
}

boolean DL_plc::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  static DL_largematrix dpdm(3,3);
  static DL_largematrix dcIdm(1,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dcIdm.resize(1,3);
    dc->dpdM(&pd,&dpdX);
    if (veloterms) {
      dc->ddpdM(&pd,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdp.times(&dpdm,&dcIdm);
    dcdm->setsubmatrix(0,0,&dcIdm);
    
    dc->dvdM(&v0,&dpdX);
    if (veloterms) {
      dc->ddvdM(&v0,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcIdm.resize(2,3);
    dcdv0.times(&dpdm,&dcIdm);
    dcdm->setsubmatrix(1,0,&dcIdm);
    
    return TRUE;
  }
  if (dc==g) { // there is an effect on C through g
    dcIdm.resize(1,3);
    dc->dpdM(&pg,&dpdX);
    if (veloterms) {
      dc->ddpdM(&pg,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdp.times(&dpdm,&dcIdm);
    dcIdm.neg();
    dcdm->setsubmatrix(0,0,&dcIdm);
    
    dc->dvdM(&w1,&dpdX);
    if (veloterms) {
      dc->ddvdM(&w1,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdw.times(&dpdm,&dcIdm);
    dcdm->setsubmatrix(1,0,&dcIdm);
    
    dc->dvdM(&w2,&dpdX);
    if (veloterms) {
      dc->ddvdM(&w2,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdw.times(&dpdm,&dcIdm);
    dcdm->setsubmatrix(2,0,&dcIdm);
    
    return TRUE;
  }
  return FALSE;
}

boolean DL_plc::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  static DL_largematrix dpdi(3,3);
  static DL_largematrix dcIdi(1,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==d) { // there is an effect on C through d
    dcIdi.resize(1,3);
    dc->dpdi(&pd,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&pd,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdp.times(&dpdi,&dcIdi);
    dcdi->setsubmatrix(0,0,&dcIdi);
    
    dc->dvdi(&v0,pc,&dpdX);
    if (veloterms) {
      dc->ddvdi(&v0,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcIdi.resize(2,3);
    dcdv0.times(&dpdi,&dcIdi);
    dcdi->setsubmatrix(1,0,&dcIdi);
    
    return TRUE;
  }
  if (dc==g) { // there is an effect on C through g
    dcIdi.resize(1,3);
    dc->dpdi(&pg,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&pg,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdp.times(&dpdi,&dcIdi);
    dcIdi.neg();
    dcdi->setsubmatrix(0,0,&dcIdi);
    
    dc->dvdi(&w1,pc,&dpdX);
    if (veloterms) {
      dc->ddvdi(&w1,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdw.times(&dpdi,&dcIdi);
    dcdi->setsubmatrix(1,0,&dcIdi);
    
    dc->dvdi(&w2,pc,&dpdX);
    if (veloterms) {
      dc->ddvdi(&w2,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdw.times(&dpdi,&dcIdi);
    dcdi->setsubmatrix(2,0,&dcIdi);
    
    return TRUE;
  }
  return FALSE;
}


void DL_plc::new_frame(void) {
 if (!g_is_dyna) {
    if (g) {
      g->new_toworld(&pg,&pgw);
      g->get_newvelocity(&pg,&dpgw);
      g->new_toworld(&w1,&w1w);
      g->get_newvelocity(&w1,&dw1w);
      g->new_toworld(&w2,&w2w);
      g->get_newvelocity(&w2,&dw2w);
    }
    else {
      pg.minus(&pgw,&dpgw);
      pgw.assign(&pg);
      w1.minus(&w1w,&dw1w);
      w1w.assign(&w1);
      w2.minus(&w2w,&dw2w);
      w2w.assign(&w2);
    }
  }

  {
    DL_vector oldl;
    oldl.assign(&l);
    d->to_world(&v0,&l);
    rotatebase(&oldl,&l,&x,&y);
  }
  
  dfdR.setcolumn(0,&l);
  dmdR.setcolumn(0,&x);
  dmdR.setcolumn(1,&y);
  dcdp.setrow(0,&l);
  dcdv0.setrow(0,&y);
  dcdv0.setrow(1,&x);
  dcdw.setrow(0,&l);

  DL_constraint::new_frame();
}

void DL_plc::test_restriction_changes(DL_largevector* lv) {
  if (maxforce>0) {
    DL_Scalar Fnew=F->get(0)+lv->get(0);
    if (Fnew>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: plane-constraint deactivated\n");
    }
  }
  if (maxtorque>0) {
    register DL_Scalar tmp1=F->get(1)+lv->get(1);
    register DL_Scalar tmp2=F->get(2)+lv->get(2);
    DL_Scalar Fnew=tmp1*tmp1+tmp2*tmp2;
    if (Fnew>maxtorque*maxtorque) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction torque: plane-constraint deactivated\n");
    }
  }
}

boolean DL_plc::check_restrictions() {
  if (maxforce>0) {
    if (F->get(0)>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: plane-constraint deactivated\n");
      return FALSE;
    }
  }
  if (maxtorque>0) {
    if ((F->get(1)*F->get(1)+F->get(2)*F->get(2))>maxtorque*maxtorque) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction torque: plane-constraint deactivated\n");
      return FALSE;
    }
  }
  return TRUE;
}

void DL_plc::apply_restrictions(DL_largevector* lv) {
  DL_vector force,torque;
  getforce(lv,&force);
  gettorque(lv,&torque);
  d->applyforce(&pd,d,&force);
  d->applytorque(&torque);
  if (g_is_dyna) {
    force.neg(&force);
    torque.neg(&torque);
    ((DL_dyna*)g)->applyforce(&pg,g,&force);
    ((DL_dyna*)g)->applytorque(&torque);
  }
}

void DL_plc::get_error(DL_largevector* lv) {
   DL_point pdw;
   DL_vector pdiff, vdiff, v0w;
   DL_Scalar hh=DL_dsystem->get_integrator()->halfstepsize();
   
   d->new_toworld(&pd,&pdw);
   if (g_is_dyna) g->new_toworld(&pg,&pgw);
   pdw.minus(&pgw,&pdiff);
   if (!testing) pdiff.timesis(stiffness);

   d->new_toworld(&v0,&v0w);
   if (g_is_dyna) {
     g->new_toworld(&w1,&w1w);
     g->new_toworld(&w2,&w2w);
   }
   if (!testing) {
     v0w.timesis(stiffness);
     w1w.timesis(stiffness);
     w2w.timesis(stiffness);
   }
   if (veloterms) {
     DL_vector dpdw,dv0w, s1, s2;

     d->get_newvelocity(&pd,&dpdw);
     if (g_is_dyna) g->get_newvelocity(&pg,&dpgw);
     dpdw.minus(&dpgw,&vdiff);
     vdiff.timesis(hh);
     pdiff.plusis(&vdiff);
     lv->set(0,pdiff.inprod(&l));
     
     d->get_newvelocity(&v0,&dv0w);
     dv0w.timesis(hh);
     if (g_is_dyna) {
       g->get_newvelocity(&w1,&dw1w);
       g->get_newvelocity(&w2,&dw2w);
       dw1w.timesis(hh); dw2w.timesis(hh);
     }
     v0w.plusis(&dv0w);
     w1w.plus(&dw1w,&s1);
     w2w.plus(&dw2w,&s2);
     lv->set(1,v0w.inprod(&s1));
     lv->set(2,v0w.inprod(&s2));
   }
   else {
     lv->set(0,pdiff.inprod(&l));
     lv->set(1,v0w.inprod(&w1w));
     lv->set(2,v0w.inprod(&w2w));
   }
}

void DL_plc::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=2;
  if (g_is_dyna) tf=4; else tf=2;
}

void DL_plc::get_force_info(int i, DL_actuator_type& at,
		            DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    at=force;
    _g=d;
    _p->assign(&pd);
    reactionforce(_v);
    return;
  }
  if (i==1) {
    at=torque;
    _g=d;
    _p->init(0,0,0);
    reactiontorque(_v);
    return;
  }
  if (g_is_dyna) {
    if (i==2) {
      at=force;
      _g=(DL_dyna*)g;
      _p->assign(&pg);
      reactionforce(_v);
      _v->neg(_v);
      return;
    }
    if (i==3) {
      at=torque;
      _g=(DL_dyna*)g;
      _p->init(0,0,0);
      reactiontorque(_v);
      _v->neg(_v);
      return;
    }
  }
  at=none;
};
