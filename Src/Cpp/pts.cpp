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
// filename     : pts.cpp
// description	: non-inline methods of class DL_pts
// author	: Bart Barenbrug     June '96
//

#include "dyna_system.h"
#include "pts.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_pts::DL_pts():DL_constraint(),
  dcdp(1,3), dfdR(3,1) {
  dim=1;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  maxforce=0.0;
  surf=NULL;
  g=NULL;
  g_is_dyna=sg_is_dyna=FALSE;
  xxinv=yyinv=0.0;
}

void DL_pts::init(DL_geo* _g, DL_point* _p, DL_surface* _surf) {
  DL_point ptmp;
  DL_vector vtmp;
  
  if (_g==_surf->get_geo()) {
    DL_dsystem->get_companion()->Msg("Error: pts::init: the surface and the point should lie in different objects!\n pts constraint not initialised\n");
    return;
  }
    
  if (_g) g_is_dyna=_g->is_dyna();
  else g_is_dyna=FALSE;
  if (_surf->get_geo()) sg_is_dyna=_surf->get_geo()->is_dyna();
  else sg_is_dyna=FALSE;
  if (!(g_is_dyna || sg_is_dyna)) {
    DL_dsystem->get_companion()->Msg("Error: pts::init: either the point or the surface should belong to a dyna\n pts-constraint not initialised\n");
    return;
  }

  // get the initial surface-parameters
  if (_g) _g->to_world(_p,&ptmp);
  else ptmp.assign(_p);
  if (!_surf->closeto(&ptmp,s,t)) {
    DL_dsystem->get_companion()->Msg("Error: pts::init: initial surface position (%f,%f) is out of bounds\n pts-constraint not initialised\n",s,t);
    return;
  }
  olds=s; oldt=t;

  // set up the local coordinate system:
  
  st_inbounds=_surf->pos(s,t,&sstl);

  if (sg_is_dyna) {
    sf=s; tf=t;
    sstf.assign(&sstl);
  }

  if (_surf->get_geo()) {
    _surf->deriv0(s,t,&vtmp);
    (_surf->get_geo())->to_world(&vtmp,&x);
  }
  else _surf->deriv0(s,t,&x);
  xxinv=1.0/x.inprod(&x);
  
  if (_surf->get_geo()) {
    _surf->deriv1(s,t,&vtmp);
    (_surf->get_geo())->to_world(&vtmp,&y);
  }
  else _surf->deriv1(s,t,&y);
  yyinv=1.0/y.inprod(&y);
  
  // make sure y is perpendicular to x:
  {
    DL_vector vtmp;
    x.times(-x.inprod(&y)*xxinv,&vtmp);
    y.plusis(&vtmp);
  }
  x.crossprod(&y,&n);
  n.normalize();

  dcdp.setrow(0,&n);
  dfdR.setcolumn(0,&n);

  surf=_surf;
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
  DL_largevector lv(1);
  get_error(&lv);
  
  if (!st_inbounds) {
    DL_dsystem->get_companion()->Msg("Error: pts::init: initial surface parameters (%f,%f) out of bounds\n pts-constraint not initialised\n",s,t);
    return;
  }
  
  if (lv.get(0)>DL_constraints->max_error) {
    DL_dsystem->get_companion()->Msg("Warning: initially invalid pts-constraint\n  Error: %f\n Initial surface parameters: (%f,%f)\n", lv.get(0), s, t);
  }

  // ok: we're in business:
  
  DL_constraint::init();
  F->makezero(); oldF->makezero();
}

void DL_pts::getforce(DL_largevector* lv, DL_vector *force) {
  n.times(lv->get(0),force);
}

void DL_pts::reactionforce(DL_vector *force) {
  getforce(F,force);
}

DL_Scalar DL_pts::get_s() {
  return s;
}

DL_Scalar DL_pts::get_t() {
  return t;
}

void DL_pts::begin_test(void) {
  DL_constraint::begin_test();
  ssave=s;
  tsave=t;
  if (sg_is_dyna)((DL_dyna *)(surf->get_geo()))->begintest();
  if (g_is_dyna) ((DL_dyna *)g)->begintest();
}

void DL_pts::end_test(void) {
  DL_constraint::end_test();
  s=ssave;
  t=tsave;
  if (sg_is_dyna)((DL_dyna *)(surf->get_geo()))->endtest();
  if (g_is_dyna) ((DL_dyna *)g)->endtest();
}


boolean DL_pts::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  static DL_largematrix dcdf;
  dcdf.resize(cc->dim,3);
  boolean nonzero=FALSE;
  if (g_is_dyna) nonzero=cc->dCdFq((DL_dyna*)g,&p,&dcdf);
  if (sg_is_dyna) {
    if (nonzero) {
      static DL_largematrix dcdfc;
      dcdfc.resize(cc->dim,3);
      if (cc->dCdFq((DL_dyna*)(surf->get_geo()),&sstl,&dcdfc))
        dcdf.minusis(&dcdfc);
    }
    else {
      if (nonzero=cc->dCdFq((DL_dyna*)(surf->get_geo()),&sstl,&dcdf))
         dcdf.neg();
    }
  }
  if (nonzero) dcdf.times(&dfdR,sub);
  return nonzero;
}

boolean DL_pts::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
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
  if (dc==surf->get_geo()) { // which implies sg_is_dyna;
               // there is an effect on C through the surface
    dc->dpdfq(&sstf,pc,&dpdX);
    if (veloterms) {
      dc->ddpdfq(&sstf,pc,&ddpdX);
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

boolean DL_pts::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
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
  if (dc==surf->get_geo()) { // which implies cg_is_dyna;
               // there is an effect on C through the surface
    dc->dpdF(&sstf, &dpdX);
    if (veloterms) {
      dc->ddpdF(&sstf, &ddpdX);
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

boolean DL_pts::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dpdm==0))
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
  if (dc==surf->get_geo()) { // which implies cg_is_dyna;
               // there is an effect on C through the curve
    dc->dpdM(&sstf, &dpdX);
    if (veloterms) {
      dc->ddpdM(&sstf, &ddpdX);
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

boolean DL_pts::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
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
  if (dc==surf->get_geo()) { // which implies sg_is_dyna;
               // there is an effect on C through the surface
    dc->dpdi(&sstf,pc,&dpdX);
    if (veloterms) {
      dc->ddpdi(&sstf,pc,&ddpdX);
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


void DL_pts::new_frame(void) {
   DL_vector vtmp;

   surf->pos(s,t,&sstl);

   if (surf->get_geo()) {
     surf->deriv0(s,t,&vtmp);
     (surf->get_geo())->to_world(&vtmp,&x);
     surf->deriv1(s,t,&vtmp);
     (surf->get_geo())->to_world(&vtmp,&y);
   }
   else {
     surf->deriv0(s,t,&x);
     surf->deriv1(s,t,&y);
   }
   xxinv=1.0/x.inprod(&x);
   yyinv=1.0/y.inprod(&y);
  
   // make sure y is perpendicular to x:
   {
     DL_vector vtmp;
     x.times(-x.inprod(&y)*xxinv,&vtmp);
     y.plusis(&vtmp);
   }
   x.crossprod(&y,&n);
   n.normalize();

   {
     DL_Scalar deltas=s-olds; olds=s; s+=deltas;
     DL_Scalar deltat=t-oldt; oldt=t; t+=deltat;
     if (sg_is_dyna) {
       DL_Scalar ast=DL_dsystem->get_integrator()->ast();
       sf=s+ast*deltas;
       tf=t+ast*deltat;
       surf->pos(sf,tf,&sstf);
     }
   }
    
   dcdp.setrow(0,&n);
   dfdR.setcolumn(0,&n);

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

void DL_pts::apply_restrictions(DL_largevector* lv) {
  DL_vector force;
  getforce(lv,&force);
  if (g_is_dyna) ((DL_dyna *)g)->applyforce(&p,g,&force);
  if (sg_is_dyna) {
    force.neg(&force);
    ((DL_dyna*)(surf->get_geo()))->applyforce(&sstf,surf->get_geo(),&force);
  }
}

void DL_pts::test_restriction_changes(DL_largevector* lv) {
  static DL_largevector Fnew(dim);
  if (!st_inbounds) {
    deactivate();
    // possibly raise an event here
  }
  if (maxforce>0) {
    F->plus(lv,&Fnew);
    if (Fnew.norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: pts-constraint deactivated\n");
    }
  }
}

boolean DL_pts::check_restrictions() {
  if (!st_inbounds) {
    deactivate();
    // possibly raise an event here
    return FALSE;
  }
  if (maxforce>0) {
    if (F->norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: pts-constraint deactivated\n");
      return FALSE;
    }
  }
  return TRUE;
}

void DL_pts::get_error(DL_largevector* lv) {
   DL_point psw, psl;
   DL_vector pdiff,vdiff,dpsw;

   // first shift s and t:
   if (surf->get_geo()) {
     st_inbounds=surf->pos(s,t,&psl);
     surf->get_geo()->new_toworld(&psl,&psw);
   }
   else st_inbounds=surf->pos(s,t,&psw);
   
   if (g_is_dyna) g->new_toworld(&p,&pgw);
    
   pgw.minus(&psw,&pdiff);
   s+=pdiff.inprod(&x)*xxinv;
   t+=pdiff.inprod(&y)*yyinv;

   // then calculate the error:
   if (!testing) pdiff.timesis(stiffness);
   if (veloterms) {
     if (surf->get_geo()) {
       psl.minus(&sstl,&vdiff);
       surf->get_geo()->new_toworld(&vdiff,&dpsw);
     }
     else psl.minus(&sstl,&dpsw);
     // velocity component due to shift in s,t is already accounted for in dpsw
     if (g_is_dyna) g->get_newvelocity(&p,&dpgw);
     dpgw.minus(&dpsw,&vdiff);
     vdiff.timesis(DL_dsystem->get_integrator()->halfstepsize());
     pdiff.plusis(&vdiff);
   }
   lv->init(pdiff.inprod(&n));
}

void DL_pts::get_fd_info(int& nrf, int& t_f) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  if (g_is_dyna) nrf=1; else nrf=0;
  if (sg_is_dyna) t_f=nrf+1; else t_f=nrf;
}

void DL_pts::get_force_info(int i, DL_actuator_type& at,
			    DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    if (g_is_dyna) {
      at=force;
      _g=(DL_dyna*)g;
      _p->assign(&p);
      reactionforce(_v);
      return;
    }
    if (sg_is_dyna) {
      at=force;
      _g=(DL_dyna*)(surf->get_geo());
      _p->assign(&sstl);
      reactionforce(_v);
      _v->neg(_v);
      return;
    }
  }
  if ((i==1) && g_is_dyna && sg_is_dyna) {
    at=force;
    _g=(DL_dyna*)(surf->get_geo());
    _p->assign(&sstl);
    reactionforce(_v);
    _v->neg(_v);
    return;
  }
  at=none;
};
