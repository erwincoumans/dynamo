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
// filename     : wheel.cpp
// description	: non-inline methods of class DL_wheel
// author	: Bart Barenbrug     October '96
//

#include "dyna_system.h"
#include "wheel.h"

// defines: WA2 causes the correction of the attachment point at the wheel
// (by linear estimate of the position halfway the next frame)
#define WA2
// CIRCLEPATH defines that the path of the wheel along the surface should
// be a circle, where the radius is determined by the rolling velocity and
// the angle between the wheel and the surface. This estimation isn't correct
// for wheels thatare attached to something else, so it is disabled
//#define CIRCLEPATH

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_wheel::DL_wheel():DL_constraint() {
  dim=3;
  F->resize(dim); F->makezero();
  oldF->resize(dim); oldF->makezero();
  Fsave->resize(dim);
  maxforce=0.0;
  w=NULL;
  surf=NULL;
  sg_is_dyna=FALSE;
  max_osc=4;
}

void DL_wheel::init(DL_dyna* _w, DL_vector* _wn, DL_point *_wc, DL_Scalar _r,
                        DL_surface *_surf, DL_Scalar _s, DL_Scalar _t) {
  DL_vector spl,vtmp,snl;
  DL_point ptmp;

  if (_w==_surf->get_geo()) {
    DL_dsystem->get_companion()->Msg("Error: wheel::init: the wheel and the surface can not be part of the same dyna\n wheel-constraint not initialised\n");
    return;
  }
  if (!_surf->indomain(_s,_t)) {
    DL_dsystem->get_companion()->Msg("Error: wheel::init: the initial surface parameters are not in the domain of the surface\n wheel-constraint not initialised\n");
    return;
  }
  st_inbounds=TRUE;
  DL_constraint::init();
  w=_w;
  surf=_surf;
  if (surf->get_geo()) sg_is_dyna=surf->get_geo()->is_dyna();
  else sg_is_dyna=FALSE;
  wr=_r;
  olds=s=_s;
  oldt=t=_t;
  w->to_local(_wn,NULL,&wn);
  wn.normalize();
  wc.assign(_wc);

  if (surf->get_geo()) {
    surf->pos(s,t,&sp);
    surf->get_geo()->to_world(&sp,&spw);
    surf->get_geo()->get_velocity(&sp,&dspw);
  }
  else {
    surf->pos(s,t,&spw);
    dspw.init(0,0,0);
  }
  
  if (surf->get_geo()) {
    surf->deriv0(s,t,&sx);
    surf->deriv1(s,t,&sy);
    sx.crossprod(&sy,&snl);
    surf->get_geo()->to_world(&snl,&snw);
    surf->get_geo()->to_world(&sx,&sxw);
    surf->get_geo()->to_world(&sy,&syw);
  }
  else {
    surf->deriv0(s,t,&sxw);
    surf->deriv1(s,t,&syw);
    sxw.crossprod(&syw,&snw);
  }
  snw.normalize();

  w->to_local(&snw,NULL,&snl);

  snl.crossprod(&wn,&wx);
  wx.normalize();
  wn.crossprod(&wx,&wy);
  wy.normalize();

  w->to_local(&spw,NULL,&ptmp);
  ptmp.minus(&wc,&spl);
  spl.normalize();

  wa=acos(spl.inprod(&wx));
  if (spl.inprod(&wy)<0) wa=-wa;
  oldwa=wa;
  DL_point wcw;
  w->to_world(&wc,&wcw);
  wcw.minus(&spw,&vtmp);
  up=(snw.inprod(&vtmp)>=0);
    
  F->makezero(); oldF->makezero();

}

void DL_wheel::reactionforce(DL_vector *f) {
  f->init(F->get(0), F->get(1), F->get(2));
}

void DL_wheel::begin_test(void) {
  DL_constraint::begin_test();
  w->begintest();
  if (sg_is_dyna) ((DL_dyna *)(surf->get_geo()))->begintest();
  ssave=s;
  tsave=t;
  wasave=wa;
}

void DL_wheel::end_test(void) {
  DL_constraint::end_test();
  w->endtest();
  if (sg_is_dyna) ((DL_dyna*)(surf->get_geo()))->endtest();
  s=ssave;
  t=tsave;
  wa=wasave;
}

boolean DL_wheel::dCdRsub(DL_constraint *cc, DL_largematrix *sub){
// return the submatrix of dCdR that shows the relation between
// the restriction value of this constraint and the constraint
// error of the constraint supplied as parameter
// the dimension of sub is cc->dim x dim;
// Returns if there is any effect at all (!result=>(m==0))
  boolean nonzero=cc->dCdFq(w,&wpc,sub);
  if (sg_is_dyna) {
    if (nonzero) {
      DL_largematrix *subtemp=new DL_largematrix(cc->dim,3);
      if (cc->dCdFq((DL_dyna*)(surf->get_geo()),&spc,subtemp))
            sub->minusis(subtemp);
      delete subtemp;
    }
    else {
      if (nonzero=cc->dCdFq((DL_dyna*)(surf->get_geo()),&spc,sub)) sub->neg();
    }
  }
  return nonzero;
}

boolean DL_wheel::dCdFq(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdfq) {
// calculate the matrix that shows the effect of application of a
// (non-central) force to the point of the dyna on the constraint
// error of this constraint.
// Returns if there is any effect at all (!result=>(dcdfq==0))
// dcdfq has dimensions dim x 3
  if (w==dc) {
    DL_matrix dpdFq, ddpdFq;
#ifdef WA2
    dc->dpdfq(&wpc2,pc,&dpdFq);
#else
    dc->dpdfq(&wpc,pc,&dpdFq);
#endif
    if (veloterms) {
#ifdef WA2
      dc->ddpdfq(&wpc2,pc,&ddpdFq);
#else
      dc->ddpdfq(&wpc,pc,&ddpdFq);
#endif
      ddpdFq.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdFq);
      dpdFq.plus(&ddpdFq,&dpdFq);
    }
    dcdfq->assign(&dpdFq);
    return TRUE;
  }
  if (surf->get_geo()==dc) {
    DL_matrix dpdFq, ddpdFq;
#ifdef WA2
    dc->dpdfq(&spc2,pc,&dpdFq);
#else
    dc->dpdfq(&spc,pc,&dpdFq);
#endif
    if (veloterms) {
#ifdef WA2
      dc->ddpdfq(&spc2,pc,&ddpdFq);
#else
      dc->ddpdfq(&spc,pc,&ddpdFq);
#endif
      ddpdFq.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdFq);
      dpdFq.plus(&ddpdFq,&dpdFq);
    }
    dcdfq->assign(&dpdFq);
    dcdfq->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_wheel::dCdF(DL_dyna *dc, DL_largematrix *dcdf) {
// calculate the matrix that shows the effect of application of a
// central force to the dyna on the constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdf==0))
// dcdf has dimensions dim x 3
  if (w==dc) {
    DL_matrix dpdf, ddpdf;
#ifdef WA2
    dc->dpdF(&wpc2,&dpdf);
#else
    dc->dpdF(&wpc,&dpdf);
#endif
    if (veloterms) {
#ifdef WA2
      dc->ddpdF(&wpc2,&ddpdf);
#else
      dc->ddpdF(&wpc,&ddpdf);
#endif
      ddpdf.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdf);
      dpdf.plus(&ddpdf,&dpdf);
    }
    dcdf->assign(&dpdf);
    return TRUE;
  }
  if (surf->get_geo()==dc) {
    DL_matrix dpdf, ddpdf;
#ifdef WA2
    dc->dpdF(&spc2,&dpdf);
#else
    dc->dpdF(&spc,&dpdf);
#endif
    if (veloterms) {
#ifdef WA2
      dc->ddpdF(&spc2,&ddpdf);
#else
      dc->ddpdF(&spc,&ddpdf);
#endif
      ddpdf.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdf);
      dpdf.plus(&ddpdf,&dpdf);
    }
    dcdf->assign(&dpdf);
    dcdf->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_wheel::dCdM(DL_dyna *dc, DL_largematrix *dcdm) {
// calculate the matrix that shows the effect of
// application of a torque to the dyna on the
// constraint error of this constraint.
// Returns if there is any effect at all (!result=>(dcdm==0))
// dcdm has dimensions dim x 3
  if (w==dc) {
    DL_matrix dpdm, ddpdm;
#ifdef WA2
    dc->dpdM(&wpc2,&dpdm);
#else
    dc->dpdM(&wpc,&dpdm);
#endif
    if (veloterms) {
#ifdef WA2
      dc->ddpdM(&wpc2,&ddpdm);
#else
      dc->ddpdM(&wpc,&ddpdm);
#endif
      ddpdm.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdm);
      dpdm.plus(&ddpdm,&dpdm);
    }
    dcdm->assign(&dpdm);
    return TRUE;
  }
  if (surf->get_geo()==dc) {
    DL_matrix dpdm, ddpdm;
#ifdef WA2
    dc->dpdM(&spc2,&dpdm);
#else
    dc->dpdM(&spc,&dpdm);
#endif
    if (veloterms) {
#ifdef WA2
      dc->ddpdM(&spc2,&ddpdm);
#else
      dc->ddpdM(&spc,&ddpdm);
#endif
      ddpdm.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdm);
      dpdm.plus(&ddpdm,&dpdm);
    }
    dcdm->assign(&dpdm);
    dcdm->neg();
    return TRUE;
  }
  return FALSE;
}

boolean DL_wheel::dCdI(DL_dyna *dc, DL_point *pc, DL_largematrix *dcdi) {
// calculate the matrix that shows the effect of application of an
// impulse to the point of the dyna on the constraint error of this
// constraint.
// Returns if there is any effect at all (!result=>(dcdi==0))
// dcdi has dimensions dim x 3
  if (w==dc) {
    DL_matrix dpdi, ddpdi;
#ifdef WA2
    dc->dpdi(&wpc2,pc,&dpdi);
#else
    dc->dpdi(&wpc,pc,&dpdi);
#endif
    if (veloterms) {
#ifdef WA2
      dc->ddpdi(&wpc2,pc,&ddpdi);
#else
      dc->ddpdi(&wpc,pc,&ddpdi);
#endif
      ddpdi.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdi);
      dpdi.plus(&ddpdi,&dpdi);
    }
    dcdi->assign(&dpdi);
    return TRUE;
  }
  if (surf->get_geo()==dc) {
    DL_matrix dpdi, ddpdi;
#ifdef WA2
    dc->dpdi(&spc2,pc,&dpdi);
#else
    dc->dpdi(&spc,pc,&dpdi);
#endif
    if (veloterms) {
#ifdef WA2
      dc->ddpdfq(&spc2,pc,&ddpdi);
#else
      dc->ddpdfq(&spc,pc,&ddpdi);
#endif
      ddpdi.times(DL_dsystem->get_integrator()->halfstepsize(),&ddpdi);
      dpdi.plus(&ddpdi,&dpdi);
    }
    dcdi->assign(&dpdi);
    dcdi->neg();
    return TRUE;
  }
  return FALSE;
}

void DL_wheel::new_frame(void) {
   DL_vector vtmp;

   wx.times(wr*cos(wa),&vtmp);
   wc.plus(&vtmp,&wpc);
   wy.times(wr*sin(wa),&vtmp);
   wpc.plus(&vtmp,&wpc);

   st_inbounds=surf->pos(s,t,&spc);

#ifdef WA2
   DL_Scalar wa2,s2,t2,ast=DL_dsystem->get_integrator()->ast();
   wa2=wa+ast*(wa-oldwa);
   s2=s+ast*(s-olds);
   t2=t+ast*(t-oldt);

   wx.times(wr*cos(wa2),&vtmp);
   wc.plus(&vtmp,&wpc2);
   wy.times(wr*sin(wa2),&vtmp);
   wpc2.plus(&vtmp,&wpc2);
   
   st_inbounds=surf->pos(s2,t2,&spc2);
#endif

   if (surf->get_geo()) {
     DL_vector snl;
     surf->deriv0(s,t,&sx);
     surf->deriv1(s,t,&sy);
     sx.crossprod(&sy,&snl);
     surf->get_geo()->to_world(&snl,&snw);
     surf->get_geo()->to_world(&sx,&sxw);
     surf->get_geo()->to_world(&sy,&syw);
   }
   else {
     surf->deriv0(s,t,&sxw);
     surf->deriv1(s,t,&syw);
     sxw.crossprod(&syw,&snw);
   }
   snw.normalize();

   // calculate rrs en rrt
   DL_vector rr;
   {
   DL_vector dwp;
   wx.times(-wr*sin(wa),&dwp);
   wy.times(wr*cos(wa),&vtmp);
   dwp.plusis(&vtmp);
   w->to_world(&dwp,&rr);
   }
   snw.times(-rr.inprod(&snw),&vtmp);
   rr.plusis(&vtmp);
   rr.normalize();

   w->to_world(&wn,&vtmp);
   DL_Scalar ftmp=vtmp.inprod(&snw);

   // if the angle between the wheel's normal and the surface normal
   // is very small (ftmp close to 1), then the wheel is almost flat
   // on the surface. This constraint type is not made for such situation,
   // so then deactivate this constraint:

   if (ftmp>0.97) {// angle smaller than ~14 degrees
     DL_dsystem->get_companion()->Msg("Wheel in wheel-constraint too flat: deactivating constraint\n");
     deactivate();
     return;
   }

//#define THICKWHEEL
#ifdef THICKWHEEL
   // apply correction for sideways elipsoid to wpc:
   DL_vector a;
   DL_Scalar beta=-ftmp;
   // limit the elipsoid angle to a ten degree maximum:
#define MAXANG 0.1
   if (beta<-MAXANG) beta=-MAXANG;
   if (beta>MAXANG) beta=MAXANG;
#undef MAXANG
#ifdef WA2
   wpc2.minus(&wc,&a); a.normalize();
   wn.times(wr*sin(beta),&vtmp);
   wpc2.plus(&vtmp,&wpc);
   a.times(wr*(cos(beta)-1),&vtmp);
   wpc2.plus(&vtmp,&wpc);
#else
   wpc.minus(&wc,&a); a.normalize();
   wn.times(wr*sin(beta),&vtmp);
   wpc.plus(&vtmp,&wpc);
   a.times(wr*(cos(beta)-1),&vtmp);
   wpc.plus(&vtmp,&wpc);
#endif
#endif

   {
   DL_point wcw,spcw;     
   w->to_world(&wc,&wcw);
   if (surf->get_geo()) surf->get_geo()->to_world(&spc,&spcw);
   else spcw.assign(&spc);
   wcw.minus(&spcw,&vtmp);
   }
   vtmp.normalize();

   DL_vector angv;
   angv.assign(w->get_angvelocity());
   angv.normalize();

#ifdef CIRCLEPATH
   straight=(fabs(angv.inprod(&vtmp))<0.1);
   
   if (straight) {
#endif
     rrs=wr*rr.inprod(&sxw)/sxw.inprod(&sxw);
     rrt=wr*rr.inprod(&syw)/syw.inprod(&syw);
#ifdef CIRCLEPATH
   }
   else {
     DL_Scalar rs;
     DL_vector vtmp2;
   
     rr.crossprod(&snw,&vtmp2);
     vtmp2.normalize();
     rs=fabs(wr*angv.inprod(&vtmp2)/angv.inprod(&vtmp));

     wrsr=wr/rs;

     sa=acos(rr.inprod(&syw)/syw.norm());
     if (rr.inprod(&sxw)>0) sa=-sa;

     rrs=rs/sxw.norm();
     rrt=rs/syw.norm();

     if ((ftmp<0)==up) {
       sa+=3.141592654;
       wrsr=-wrsr;
     }
   }
#endif
      
   DL_Scalar delta=wa-oldwa; oldwa=wa; wa+=delta;
   olds=s; oldt=t;
#ifdef CIRCLEPATH
   if (straight) {
#endif
      s+=rrs*delta; t+=rrt*delta;
#ifdef CIRCLEPATH
   }
   else {
      DL_Scalar saold=sa;
      sa+=delta*wrsr;
      s+=rrs*(cos(sa)-cos(saold));
      t+=rrt*(sin(sa)-sin(saold));
   }
#endif

   DL_constraint::new_frame();
}

void DL_wheel::apply_restrictions(DL_largevector* lv) {
  DL_vector force;
  force.init(lv->get(0), lv->get(1), lv->get(2));
#ifdef WA2
  w->applyforce(&wpc2,w,&force);
#else
  w->applyforce(&wpc,w,&force);
#endif
  if (sg_is_dyna) {
    force.neg(&force);
#ifdef WA2
    ((DL_dyna*)(surf->get_geo()))->applyforce(&spc2,surf->get_geo(),&force);
#else
    ((DL_dyna*)(surf->get_geo()))->applyforce(&spc,surf->get_geo(),&force);
#endif
  }
}

void DL_wheel::test_restriction_changes(DL_largevector* lv) {
  if (!st_inbounds) {
    deactivate();
    // possibly raise an event here
  }
  if (maxforce>0) {
    DL_largevector* Fnew=new DL_largevector(dim);
    F->plus(lv,Fnew);
    if (Fnew->norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: wheel-constraint deactivated\n");
    }
    delete Fnew;
  }
}

boolean DL_wheel::check_restrictions() {
  if (!st_inbounds) {
    deactivate();
    // possibly raise an event here
    return FALSE;
  }
  if (maxforce>0) {
    if (F->norm()>maxforce) {
      deactivate();
      // possibly raise an event here
      DL_dsystem->get_companion()->Msg("Too large a reaction force: wheel-constraint deactivated\n");
      return FALSE;
    }
  }
  return TRUE;
}

void DL_wheel::get_error(DL_largevector* lv) {
   DL_vector dwpw, pdiff, vdiff, vtmp;
   DL_point dwp;

   wx.times(wr*cos(wa),&vdiff);
   wc.plus(&vdiff,&wp);
   wy.times(wr*sin(wa),&vdiff);
   wp.plus(&vdiff,&wp);
   w->new_toworld(&wp,&wpw);

   st_inbounds=surf->pos(s,t,&sp);
   if (surf->get_geo()) surf->get_geo()->new_toworld(&sp,&spw);
   else spw.assign(&sp);
   wpw.minus(&spw,&pdiff);
   if (!testing) pdiff.timesis(stiffness);
   if (veloterms) {
     w->get_newvelocity(&wp,&dwpw);
     if (surf->get_geo()) surf->get_geo()->get_newvelocity(&sp,&dspw);
     else dspw.init(0,0,0);
     // there's also a velocity component due to the changes in s, t and w,
     // but these components are (by construction) the same for the point
     // in the surface and in the wheel, so they cancel each other and need
     // not explicitly be accounted for.
     dwpw.minus(&dspw,&vdiff);
     vdiff.timesis(DL_dsystem->get_integrator()->halfstepsize());
     // terms for the motions of the contact points in the wheel and the
     // surface are not required since both displacements are calculated
     // from the change in wa, and therefore cancel each other
     pdiff.plusis(&vdiff);
   }
   lv->init(pdiff.x, pdiff.y, pdiff.z);

   // next adjust wa, s and t according to tangent error

   // get the tangent (in wc):
   wx.times(-wr*sin(wa),&vtmp);
   wy.times(wr*cos(wa),&vdiff);
   vtmp.plusis(&vdiff);
   w->new_toworld(&vtmp,&dwpw);
   dwpw.normalize();

   // get the surface normal (in wc):
   if (sg_is_dyna) {
     DL_vector ds,dt,snl;
     surf->deriv0(s,t,&ds);
     surf->deriv1(s,t,&dt);
     ds.crossprod(&dt,&snl);
     surf->get_geo()->new_toworld(&snl,&snw);
     snw.normalize();
   }

   DL_Scalar da=dwpw.inprod(&snw);
   if (up) da=-da;
   wa+=da;

   if (straight) {
      s+=rrs*da; t+=rrt*da;
   }
   else {
      DL_Scalar saold=sa;
      sa+=da*wrsr;
      s+=rrs*(cos(sa)-cos(saold));
      t+=rrt*(sin(sa)-sin(saold));
   }
}

void DL_wheel::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=1;
  if (sg_is_dyna) tf=2; else tf=1;
}

void DL_wheel::get_force_info(int i, DL_actuator_type& at,
			      DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    at=force;
    _g=w;
    _p->assign(&wpc);
    reactionforce(_v);
    return;
  }
  if ((i==1) && (sg_is_dyna)) {
    at=force;
    _g=(DL_dyna*)surf->get_geo();
    _p->assign(&spc);
    reactionforce(_v);
    _v->neg(_v);
    return;
  }
  at=none;
};
