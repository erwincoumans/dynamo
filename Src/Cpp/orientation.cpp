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
// filename     : orientation.cpp
// description	: non-inline methods of class orientation
// author	: Bart Barenbrug     August 1996
//

#include "orientation.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_orientation::DL_orientation():DL_constraint(),
  dcdv0(3,3),dcdv1(1,3),dcdw1(1,3),dcdw2(3,3),dfdr(3,3) {
  dim=3;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  maxtorque=0.0;
  d=NULL;
  g=NULL;
}

DL_orientation::~DL_orientation() {
}

void DL_orientation::gettorque(DL_largevector* lv,DL_vector *t) {
  DL_vector tmp;
  x.times(lv->get(0),t);
  y.times(lv->get(1),&tmp); t->plusis(&tmp);
  z.times(lv->get(2),&tmp); t->plusis(&tmp);
}

void DL_orientation::reactiontorque(DL_vector *t) {
  gettorque(F,t);
}

void DL_orientation::get_dyna_vector0(DL_vector *v){
  v->assign(&v0);
}

void DL_orientation::set_dyna_vector0(DL_vector *v){
  DL_Scalar l=v->norm();
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Warning: orientation::set_dyna_vector0(): nul-vector assignment ignored\n");
    return;
  }
  else v->times(d->torquefactor()/l,&v0);
}

void DL_orientation::get_dyna_vector1(DL_vector *v){
  v->assign(&v1);
}

void DL_orientation::set_dyna_vector1(DL_vector *v){
  DL_Scalar l=v->norm();
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Warning: orientation::set_dyna_vector1(): zero-vector assignment ignored\n");
    return;
  }
  else v->times(d->torquefactor()/l,&v1);
}

void DL_orientation::get_geo_vector1(DL_vector *v){
   v->assign(&w1);
}

void DL_orientation::set_geo_vector1(DL_vector *v){
  DL_Scalar l=v->norm();
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Warning: orientation::set_geo_vector1(): zero-vector assignment ignored\n");
    return;
  }
  else v->times(d->torquefactor()/l,&w1);
}

void DL_orientation::get_geo_vector2(DL_vector *v){
  v->assign(&w2);
}

void DL_orientation::set_geo_vector2(DL_vector *v){
  DL_Scalar l=v->norm();
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Warning: orientation::set_geo_vector2(): nul-vector assignment ignored\n");
    return;
  }
  else v->times(d->torquefactor()/l,&w2);
}

void DL_orientation::begin_test(void) {
  DL_constraint::begin_test();
  d->begintest();
  if (g_is_dyna) ((DL_dyna *)g)->begintest();
}

void DL_orientation::end_test(void) {
  DL_constraint::end_test();
  d->endtest();
  if (g_is_dyna) ((DL_dyna*)g)->endtest();
}

void DL_orientation::init(DL_dyna* _d, DL_vector *_v0, DL_vector *_v1,
                              DL_geo* _g, DL_vector *_w1, DL_vector *_w2) {
  if (_d==_g) {
    DL_dsystem->get_companion()->Msg("Error: orientation::init: an orientation constraint needs vectors from _different_ objects\n orientation constraint not initialised\n");
    return;
  }
  d=_d;
  v0.assign(_v0);
  v1.assign(_v1);
  g=_g;
  w1.assign(_w1);
  w2.assign(_w2);
  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;
  DL_Scalar tf=d->torquefactor();
  DL_Scalar l=v0.norm();
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Error: orientation constraint: direction vector v0 should be non-zero!\n orientation constraint not initialised\n");
    return;
  }
  else v0.timesis(tf/l);
  l=v1.norm();
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Error: orientation constraint: direction vector v1 should be non-zero!\n orientation constraint not initialised\n");
    return;
  }
  else v1.timesis(tf/l);
  l=w1.norm();
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Error: orientation constraint: direction vector w1 should be non-zero!\n orientation constraint not initialised\n");
    return;
  }
  else w1.timesis(tf/l);
  l=w2.norm();
  if (l==0) {
    DL_dsystem->get_companion()->Msg("Error: orientation constraint: direction vector w2 should be non-zero!\n orientation constraint not initialised\n");
    return;
  }
  else w2.timesis(tf/l);
  DL_constraint::init();
  if (!g_is_dyna) {
    if (g) {
      g->new_toworld(&w1,&w1w);
      g->new_toworld(&w2,&w2w);
      g->get_newvelocity(&w1,&dw1w);
      g->get_newvelocity(&w2,&dw2w);
    }
    else {
      w1w.assign(&w1);
      w2w.assign(&w2);
      dw1w.init(0,0,0);
      dw2w.init(0,0,0);
    }
  }
  
  F->makezero(); oldF->makezero();

  // check if the constraint is initially valid
  {
    DL_largevector lv(3);
    get_error(&lv);
    if (lv.norm()!=0.0) {
      DL_dsystem->get_companion()->Msg("Warning: initially invalid orientation constraint\n  Error: (%f,%f,%f)\n",
				    lv.get(0), lv.get(1), lv.get(2) );
    }
  }
  
  DL_vector vtmp(0,0,0);
  dcdv0.setrow(2,&vtmp);
  dcdw2.setrow(0,&vtmp);
  d->to_world(&v0,&z);
  dcdw1.setrow(0,&z);
  dcdw2.setrow(1,&z);
  z.normalize();
  dfdr.setcolumn(2,&z);
  g->to_world(&w1,&y);
  dcdv0.setrow(0,&y);
  dcdw2.setrow(2,&y);
  y.normalize();
  dfdr.setcolumn(1,&y);
  g->to_world(&w2,&x);
  dcdv0.setrow(1,&x);
  dcdv1.setrow(0,&x);
  x.normalize();
  dfdr.setcolumn(0,&x);
}

boolean DL_orientation::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  static DL_largematrix dcdf,dcdfg;
  dcdf.resize(cc->dim,3);
  boolean nonzero=cc->dCdM(d,&dcdf);
  if (g_is_dyna) {
    if (nonzero) {
      dcdfg.resize(cc->dim,3);
      if (cc->dCdM((DL_dyna*)g,&dcdfg)) dcdf.minusis(&dcdfg);
    }
    else if (nonzero=cc->dCdM((DL_dyna*)g,&dcdf)) dcdf.neg(); 
  }
  dcdf.times(&dfdr,sub);
  return nonzero;
}

boolean DL_orientation::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  static DL_largematrix dvdf(3,3);
  static DL_largematrix dcdftmp(1,3);
  DL_matrix dvdX;
  DL_matrix ddvdX;
  if (d==dc) {
    dcdfq->makezero();
    
    dc->dvdfq(&v1,pc,&dvdX);
    if (veloterms) {
      dc->ddvdfq(&v1,pc,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdf.assign(&dvdX);
    dcdftmp.resize(1,3);
    dcdv1.times(&dvdf,&dcdftmp);
    dcdfq->setsubmatrix(2,0,&dcdftmp);

    dc->dvdfq(&v0,pc,&dvdX);
    if (veloterms) {
      dc->ddvdfq(&v0,pc,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdf.assign(&dvdX);
    dcdftmp.resize(3,3);
    dcdv0.times(&dvdf,&dcdftmp);
    dcdfq->plusis(&dcdftmp);
    
    return TRUE;
  }
  if (g==dc) {
    dcdfq->makezero();
    
    dc->dvdfq(&w1,pc,&dvdX);
    if (veloterms) {
      dc->ddvdfq(&w1,pc,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdf.assign(&dvdX);
    dcdftmp.resize(1,3);
    dcdw1.times(&dvdf,&dcdftmp);
    dcdfq->setsubmatrix(0,0,&dcdftmp);

    dc->dvdfq(&w2,pc,&dvdX);
    if (veloterms) {
      dc->ddvdfq(&w2,pc,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdf.assign(&dvdX);
    dcdftmp.resize(3,3);
    dcdw2.times(&dvdf,&dcdftmp);
    dcdfq->plusis(&dcdftmp);
    
    return TRUE;
  }
  return FALSE;
}

boolean DL_orientation::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  // a central force has no effect on the orientation:	      
  return FALSE;
}

boolean DL_orientation::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  static DL_largematrix dvdm(3,3);
  static DL_largematrix dcdmtmp(1,3);
  DL_matrix dvdX;
  DL_matrix ddvdX;
  if (d==dc) {
    dcdm->makezero();
    
    dc->dvdM(&v1,&dvdX);
    if (veloterms) {
      dc->ddvdM(&v1,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdm.assign(&dvdX);
    dcdmtmp.resize(1,3);
    dcdv1.times(&dvdm,&dcdmtmp);
    dcdm->setsubmatrix(2,0,&dcdmtmp);

    dc->dvdM(&v0,&dvdX);
    if (veloterms) {
      dc->ddvdM(&v0,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdm.assign(&dvdX);
    dcdmtmp.resize(3,3);
    dcdv0.times(&dvdm,&dcdmtmp);
    dcdm->plusis(&dcdmtmp);
    
    return TRUE;
  }
  if (g==dc) {
    dcdm->makezero();
    
    dc->dvdM(&w1,&dvdX);
    if (veloterms) {
      dc->ddvdM(&w1,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdm.assign(&dvdX);
    dcdmtmp.resize(1,3);
    dcdw1.times(&dvdm,&dcdmtmp);
    dcdm->setsubmatrix(0,0,&dcdmtmp);

    dc->dvdM(&w2,&dvdX);
    if (veloterms) {
      dc->ddvdM(&w2,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdm.assign(&dvdX);
    dcdmtmp.resize(3,3);
    dcdw2.times(&dvdm,&dcdmtmp);
    dcdm->plusis(&dcdmtmp);
    
    return TRUE;
  }
  return FALSE;
}

boolean DL_orientation::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  static DL_largematrix dvdi(3,3);
  static DL_largematrix dcditmp(1,3);
  DL_matrix dvdX;
  DL_matrix ddvdX;
  if (d==dc) {
    dcdi->makezero();
    
    dc->dvdi(&v1,pc,&dvdX);
    if (veloterms) {
      DL_matrix ddvdX;
      dc->ddvdi(&v1,pc,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdi.assign(&dvdX);
    dcditmp.resize(1,3);
    dcdv1.times(&dvdi,&dcditmp);
    dcdi->setsubmatrix(2,0,&dcditmp);

    dc->dvdi(&v0,pc,&dvdX);
    if (veloterms) {
      dc->ddvdi(&v0,pc,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdi.assign(&dvdX);
    dcditmp.resize(3,3);
    dcdv0.times(&dvdi,&dcditmp);
    dcdi->plusis(&dcditmp);
    
    return TRUE;
  }
  if (g==dc) {
    dcdi->makezero();
    
    dc->dvdi(&w1,pc,&dvdX);
    if (veloterms) {
      dc->ddvdi(&w1,pc,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdi.assign(&dvdX);
    dcditmp.resize(1,3);
    dcdw1.times(&dvdi,&dcditmp);
    dcdi->setsubmatrix(0,0,&dcditmp);

    dc->dvdi(&w2,pc,&dvdX);
    if (veloterms) {
      dc->ddvdi(&w2,pc,&ddvdX);
      ddvdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddvdX);
      dvdX.plus(&ddvdX,&dvdX);
    }
    dvdi.assign(&dvdX);
    dcditmp.resize(3,3);
    dcdw2.times(&dvdi,&dcditmp);
    dcdi->plusis(&dcditmp);
    
    return TRUE;
  }
  return FALSE;
}


void DL_orientation::new_frame(void) {
  d->to_world(&v0,&z);
  dcdw1.setrow(0,&z);
  dcdw2.setrow(1,&z);
  z.normalize();
  dfdr.setcolumn(2,&z);
  d->to_world(&v1,&y);
  dcdv0.setrow(0,&y);
  dcdw2.setrow(2,&y);
  y.normalize();
  dfdr.setcolumn(1,&y);
  g->to_world(&w2,&x);
  dcdv0.setrow(1,&x);
  dcdv1.setrow(0,&x);
  x.normalize();
  dfdr.setcolumn(0,&x);
  if (!g_is_dyna) {
    if (g) {
      g->new_toworld(&w1,&w1w);
      g->new_toworld(&w2,&w2w);
      g->get_newvelocity(&w1,&dw1w);
      g->get_newvelocity(&w2,&dw2w);
    }
    else {
      DL_Scalar hinv=1.0/DL_dsystem->get_integrator()->old_stepsize();
      w1.minus(&w1w,&dw1w);
      w2.minus(&w2w,&dw2w);
      dw1w.timesis(hinv);
      dw2w.timesis(hinv);
      w1w.assign(&w1);
      w2w.assign(&w2);
    }
  }
  DL_constraint::new_frame();
}

void DL_orientation::apply_restrictions(DL_largevector* lv) {
  DL_vector torque;
  gettorque(lv,&torque);
  d->applytorque(&torque);
  if (g_is_dyna) {
    torque.neg(&torque);
    ((DL_dyna*)g)->applytorque(&torque);
  }
}

boolean DL_orientation::check_restrictions() {
  if (maxtorque>0) {
    DL_vector trq;
    gettorque(F,&trq);
    if (trq.norm()>maxtorque) {
      DL_dsystem->get_companion()->Msg("Too large a reaction torque: orientation-constraint deactivated\n");
      // possibly raise an event here
      deactivate();
      return FALSE;
    }
  }
  return TRUE;
}

void DL_orientation::test_restriction_changes(DL_largevector* lv) {
  if (maxtorque>0) {
    DL_largevector* Fnew=new DL_largevector(dim);
    F->plus(lv,Fnew);
    DL_vector trq;
    gettorque(Fnew,&trq);
    if (trq.norm()>maxtorque) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction torque: orientation-constraint deactivated\n");
    }
    delete Fnew;
  }
}

void DL_orientation::get_error(DL_largevector* lv) {
   DL_vector v0w, v1w;
   DL_Scalar hh=DL_dsystem->get_integrator()->halfstepsize();
   
   d->new_toworld(&v0,&v0w);
   d->new_toworld(&v1,&v1w);
   if (g_is_dyna) {
     g->new_toworld(&w1,&w1w);
     g->new_toworld(&w2,&w2w);
   }
   if (!testing) {
     v0w.timesis(stiffness); v1w.timesis(stiffness);
     w1w.timesis(stiffness); w2w.timesis(stiffness);
   }
   
   if (veloterms) {
     DL_vector dv0w, dv1w, s1,s2;
     d->get_newvelocity(&v0,&dv0w);
     d->get_newvelocity(&v1,&dv1w);
     if (g_is_dyna) {
       g->get_newvelocity(&w1,&dw1w);
       g->get_newvelocity(&w2,&dw2w);
       dw1w.timesis(hh); dw2w.timesis(hh);
     }   
     dv0w.timesis(hh); dv1w.timesis(hh);
     v0w.plusis(&dv0w); v1w.plusis(&dv1w);
     w1w.plus(&dw1w,&s1);
     w2w.plus(&dw2w,&s2);
     lv->init(v0w.inprod(&s1), v0w.inprod(&s2), v1w.inprod(&s2));
   }
   else lv->init(v0w.inprod(&w1w),v0w.inprod(&w2w),v1w.inprod(&w2w));
 }

void DL_orientation::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=1;
  if (g_is_dyna) tf=2; else tf=1;
}

void DL_orientation::get_force_info(int i, DL_actuator_type& at,
			      DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    at=torque;
    _g=d;
    _p->init(0,0,0);
    reactiontorque(_v);
    return;
  }
  if ((i==1) && (g_is_dyna)) {
    at=torque;
    _g=(DL_dyna*)g;
    _p->init(0,0,0);
    reactiontorque(_v);
    _v->neg(_v);
    return;
  }
  at=none;
};
