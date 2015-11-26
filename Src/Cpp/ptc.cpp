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
// filename     : ptc.cpp
// description	: non-inline methods of class DL_ptc
// author	: Bart Barenbrug     June '96
//

#include "dyna_system.h"
#include "ptc.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_ptc::DL_ptc():DL_constraint(),
  dcdp(2,3), dfdR(3,2) {
  dim=2;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  maxforce=0.0;
  c=NULL;
  g=NULL;
  g_is_dyna=cg_is_dyna=FALSE;
  ddinv=0.0;
}

void DL_ptc::init(DL_geo* _g, DL_point* _p, DL_curve* _c) {
  DL_point ptmp;

  if (_g==_c->get_geo()) {
    DL_dsystem->get_companion()->Msg("Error: ptc::init: the curve and the point should lie in different objects!\n ptc constraint not initialised\n");
    return;
  }
    
  if (_g) g_is_dyna=_g->is_dyna(); else g_is_dyna=FALSE;
  if (_c->get_geo()) cg_is_dyna=_c->get_geo()->is_dyna(); else cg_is_dyna=FALSE;
  if (!(g_is_dyna || cg_is_dyna)) {
    DL_dsystem->get_companion()->Msg("Error: ptc::init: either the point or the curve should belong to a dyna\n ptc-constraint not initialised\n");
    return;
  }

  // get the initial curve-parameter
  if (_g) {
    _g->to_world(_p,&ptmp);
    sf=olds=s=_c->closeto(&ptmp);
  }
  else olds=s=_c->closeto(_p);

  // set up the local coordinate system:
  _c->pos(s,&csl);
  if (cg_is_dyna) csf.assign(&csl);

  if (_c->get_geo()) {
    DL_vector vtmp;
    _c->deriv(s,&vtmp);
    _c->get_geo()->to_world(&vtmp,&d);
  }
  else _c->deriv(s,&d);
  ddinv=1.0/d.inprod(&d);

  // x and y are two vectors each with length 1 and both perpendicular to
  // each other and d
  if ((d.x!=0) || (d.z!=0)) x.init(d.z,0,-(d.x));
  else x.init(0,0,-(d.y));
  x.normalize();
  d.crossprod(&x,&y);
  y.normalize();

  dcdp.setrow(0,&x);
  dcdp.setrow(1,&y);
  dfdR.setcolumn(0,&x);
  dfdR.setcolumn(1,&y);

  c=_c;
  g=_g;
  p.assign(_p);

  if (!g_is_dyna) {
    if (g) {
      g->new_toworld(&p,&pgw);
      g->get_newvelocity(&p,&dpgw);
    }
    else {
      pgw.assign(&p);
      dpgw.init(0,0,0);
    }
  }

  // check if the constraint is initially valid
  DL_largevector lv(2);
  get_error(&lv);
  
  if (!s_inbounds) {
    DL_dsystem->get_companion()->Msg("Error: ptc::init: initial curve position (%f) is out of bounds\n ptc-constraint not initialised", s);
    return;
  }
  
  if (lv.norm()>DL_constraints->max_error) {
    DL_dsystem->get_companion()->Msg("Warning: initially invalid ptc-constraint\n Error: (%f,%f)\n Initial curve parameter: %s\n",
				  lv.get(0), lv.get(1), s );
  }

  // ok: we're in business:
  
  DL_constraint::init();
  F->makezero(); oldF->makezero();
}

void DL_ptc::getforce(DL_largevector* lv, DL_vector *force) {
  DL_vector vtmp;
  x.times(lv->get(0),force);
  y.times(lv->get(1),&vtmp);
  force->plusis(&vtmp);
}

void DL_ptc::reactionforce(DL_vector *force) {
  getforce(F,force);
}

DL_Scalar DL_ptc::get_s() {
  return s;
}

void DL_ptc::begin_test(void) {
  DL_constraint::begin_test();
  ssave=s;
  if (cg_is_dyna)((DL_dyna *)(c->get_geo()))->begintest();
  if (g_is_dyna) ((DL_dyna *)g)->begintest();
}

void DL_ptc::end_test(void) {
  DL_constraint::end_test();
  s=ssave;
  if (cg_is_dyna)((DL_dyna *)(c->get_geo()))->endtest();
  if (g_is_dyna) ((DL_dyna *)g)->endtest();
}

boolean DL_ptc::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  static DL_largematrix dcdf;
  static DL_largematrix dcdfc;
  boolean nonzero=FALSE;
  dcdf.resize(cc->dim,3);
  if (g_is_dyna) nonzero=cc->dCdFq((DL_dyna*)g,&p,&dcdf);
  if (cg_is_dyna) {
    if (nonzero) {
      dcdfc.resize(cc->dim,3);
      if (cc->dCdFq((DL_dyna*)(c->get_geo()),&csl,&dcdfc))
        dcdf.minusis(&dcdfc);
    }
    else {
      if (nonzero=cc->dCdFq((DL_dyna*)(c->get_geo()),&csl,&dcdf)) dcdf.neg();
    }
  }
  if (nonzero) dcdf.times(&dfdR,sub);
  return nonzero;
}

boolean DL_ptc::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  static DL_largematrix dpdfq(3,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==g) { // there is an effect on C through g
    dc->dpdfq(&p,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&p,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdfq.assign(&dpdX);
    dcdp.times(&dpdfq,dcdfq);
    return TRUE;
  }
  if (dc==c->get_geo()) { // which implies cg_is_dyna;
               // there is an effect on C through the curve
    dc->dpdfq(&csf,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&csf,pc,&ddpdX);
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

boolean DL_ptc::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  static DL_largematrix dpdi(3,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==g) { // there is an effect on C through g
    dc->dpdi(&p,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&p,pc,&ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdi.assign(&dpdX);
    dcdp.times(&dpdi,dcdi);
    return TRUE;
  }
  if (dc==c->get_geo()) { // which implies cg_is_dyna;
               // there is an effect on C through the curve
    dc->dpdi(&csf,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&csf,pc,&ddpdX);
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

boolean DL_ptc::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  static DL_largematrix dpdf(3,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==g) { // there is an effect on C through g
    dc->dpdF(&p, &dpdX);
    if (veloterms) {
      dc->ddpdF(&p, &ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdf.assign(&dpdX);
    dcdp.times(&dpdf,dcdf);
    return TRUE;
  }
  if (dc==c->get_geo()) { // which implies cg_is_dyna;
               // there is an effect on C through the curve
    dc->dpdF(&csf, &dpdX);
    if (veloterms) {
      dc->ddpdF(&csf, &ddpdX);
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

boolean DL_ptc::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  static DL_largematrix dpdm(3,3);
  DL_matrix dpdX;
  DL_matrix ddpdX;
  if (dc==g) { // there is an effect on C through g
    dc->dpdM(&p, &dpdX);
    if (veloterms) {
      dc->ddpdM(&p, &ddpdX);
      ddpdX.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdX);
      dpdX.plus(&ddpdX,&dpdX);
    }
    dpdm.assign(&dpdX);
    dcdp.times(&dpdm,dcdm);
    return TRUE;
  }
  if (dc==c->get_geo()) { // which implies cg_is_dyna;
               // there is an effect on C through the curve
    dc->dpdM(&csf, &dpdX);
    if (veloterms) {
      dc->ddpdM(&csf, &ddpdX);
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

void DL_ptc::new_frame(void) {
   DL_vector oldd,vtmp;

   c->pos(s,&csl);
   oldd.assign(&d);

   if (c->get_geo()) {
     c->deriv(s,&vtmp);
     (c->get_geo())->to_world(&vtmp,&d);
   }
   else c->deriv(s,&d);
   ddinv=1.0/d.inprod(&d);

   rotatebase(&oldd,&d,&x,&y);

   DL_Scalar deltas=s-olds;
   olds=s;
   s+=deltas;

   if (!c->indomain(s)) {
     reset();
     deactivate();
     // possibly raise an event here
     return;
   }
   
   if (cg_is_dyna) {
     DL_Scalar ast=DL_dsystem->get_integrator()->ast();
     sf=s+ast*deltas;
     c->pos(sf,&csf);
   }
     
   dcdp.setrow(0,&x);
   dcdp.setrow(1,&y);
   dfdR.setcolumn(0,&x);
   dfdR.setcolumn(1,&y);

   if (!g_is_dyna) {
     if (g) {
       g->new_toworld(&p,&pgw);
       g->get_newvelocity(&p,&dpgw);
     }
     else {
       p.minus(&pgw,&dpgw);
       dpgw.timesis(1.0/DL_dsystem->get_integrator()->old_stepsize());
       pgw.assign(&p);
     }
   }

   DL_constraint::new_frame();
}

void DL_ptc::apply_restrictions(DL_largevector* lv) {
  DL_vector force;
  getforce(lv,&force);
  if (g_is_dyna) ((DL_dyna *)g)->applyforce(&p,g,&force);
  if (cg_is_dyna) {
    force.neg(&force);
    ((DL_dyna*)(c->get_geo()))->applyforce(&csf,c->get_geo(),&force);
  }
}

void DL_ptc::get_error(DL_largevector* lv) {
   DL_point pcl, pcw;
   DL_vector pdiff, dpcw, vdiff;

   // first shift s:
   if (c->get_geo()) {
     s_inbounds=c->pos(s,&pcl);
     c->get_geo()->new_toworld(&pcl,&pcw);
   }
   else s_inbounds=c->pos(s,&pcw);
   if (g_is_dyna) g->new_toworld(&p,&pgw);
    
   pgw.minus(&pcw,&pdiff);
   s+=pdiff.inprod(&d)*ddinv;

   // then calculate the error
   if (!testing) pdiff.timesis(stiffness);
   if (veloterms) {
     if (c->get_geo()) {
       pcl.minus(&csl,&vdiff);
       c->get_geo()->new_toworld(&vdiff,&dpcw);
     }
     else pcl.minus(&csl,&dpcw);
     // velocity component due to shifting s is already incorporated in dpcw
     if (g_is_dyna) g->get_newvelocity(&p,&dpgw);
     dpgw.minus(&dpcw,&vdiff);
     vdiff.timesis(DL_dsystem->get_integrator()->halfstepsize());
     pdiff.plusis(&vdiff);
   }
   lv->init(pdiff.inprod(&x),pdiff.inprod(&y));
}

void DL_ptc::test_restriction_changes(DL_largevector* lv) {
  static DL_largevector Fnew(dim);
  if (!s_inbounds) {
    deactivate();
    // possibly raise an event here
  }
  if (maxforce>0) {
    F->plus(lv,&Fnew);
    if (Fnew.norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: ptc-constraint deactivated\n");
    }
  }
}

boolean DL_ptc::check_restrictions() {
  if (!s_inbounds) {
    deactivate();
    // possibly raise an event here
    return FALSE;
  }
  if (maxforce>0) {
    if (F->norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: ptc-constraint deactivated\n");
      return FALSE;
    }
  }
  return TRUE;
}

void DL_ptc::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  if (g_is_dyna) nrf=1; else nrf=0;
  if (cg_is_dyna) tf=nrf+1; else tf=nrf;
}

void DL_ptc::get_force_info(int i, DL_actuator_type& at,
			      DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    if (g_is_dyna) {
      at=force;
      _g=(DL_dyna*)g;
      _p->assign(&p);
      reactionforce(_v);
      return;
    }
    if (cg_is_dyna) {
      at=force;
      _g=(DL_dyna*)(c->get_geo());
      _p->assign(&csl);
      reactionforce(_v);
      _v->neg(_v);
      return;
    }
  }
  if ((i==1) && g_is_dyna && cg_is_dyna) {
    at=force;
    _g=(DL_dyna*)(c->get_geo());
    _p->assign(&csl);
    reactionforce(_v);
    _v->neg(_v);
    return;
  }
  at=none;
};
