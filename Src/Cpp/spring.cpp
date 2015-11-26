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
// filename     : spring.cpp
// description	: methods of class DL_spring
// author	: Bart Barenbrug     August '96
//

#include "dyna_system.h"
#include "spring.h"

// *************** //
// member fuctions //
// *************** //

DL_spring::DL_spring() {
  c=dc=l=maxforce=0; g_is_dyna=el=FALSE;
  d=NULL; g=NULL;
}

DL_spring::~DL_spring(void) {
}

void DL_spring::init(DL_dyna *_d, DL_point *_pd,
                            DL_geo *_g, DL_point *_pg) {
  d=_d; g=_g;
  pd.assign(_pd);
  pg.assign(_pg);
  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;
  sforce.init(0,0,0);
  activate();
}

void DL_spring::init(DL_dyna *_d, DL_point *_pd,
                            DL_geo *_g, DL_point *_pg,
			    DL_Scalar _l, DL_Scalar _c) {
  c=_c;
  l=_l;
  init(_d,_pd,_g,_pg);
}

void DL_spring::init(DL_dyna *_d, DL_point *_pd,
                            DL_geo *_g, DL_point *_pg,
			    DL_Scalar _l, DL_Scalar _c, DL_Scalar _dc) {
  dc=_dc;
  init(_d,_pd,_g,_pg,_l,_c);
}

#define DL_MIN_NORMALIZE_LEN   0.00001

boolean DL_spring::calculate_force(void) {
  DL_point pdw,pgw;
  DL_vector dpdw,dpgw;
  DL_vector pdiff,vdiff;
  DL_Scalar len,vterm,df,fac;
  d->to_world(&pd,&pdw);
  if (g) {
    g->to_world(&pg,&pgw);
    pgw.minus(&pdw,&pdiff);
  }
  else pg.minus(&pdw,&pdiff);
  len=pdiff.norm();
  if (len > DL_MIN_NORMALIZE_LEN) { // safer code added by Gary R. Van Sickle
       // pdiff is long enough to be normalized
       pdiff.normalize();
  } 
  else {
       // Not long enough, but len is so small we don't really care, so just pick a random one
       pdiff.init(1.0, 0.0, 0.0);
  }
  // account for discretisation errors (since the reaction force
  // shouldn't be constant over the next frame, but should vary
  // since the length also varies according to the relative speed
  // in the direction of spring:
  d->get_velocity(&pd,&dpdw);
  if (g) {
    g->get_velocity(&pg,&dpgw);
    dpdw.minus(&dpgw,&vdiff);
  }
  else vdiff.assign(&dpdw);
  df=0.5*DL_dsystem->get_integrator()->stepsize();  // discretisation factor
  // calculating the damping term:
  vterm=vdiff.inprod(&pdiff);
  // and combine them all into the reaction force:
  fac=c*(len-df*vterm-l)-dc*vterm;
  if (el) {
    if (fac<0) {
      sforce.init(0,0,0);
      return TRUE;
    }
  }
  pdiff.times(fac,&sforce);
  return FALSE;
}

#undef DL_MIN_NORMALIZE_LEN

void DL_spring::calculate_and_apply(void) {
  if (calculate_force()) return;
  if (maxforce!=0) {
    if (maxforce>0) {
      if (sforce.norm()>maxforce) {
        deactivate();
        // possibly raise an event here
        return;
      }
    }
    else { // maxforce<0
      DL_Scalar fl=-sforce.norm();
      if (fl<maxforce) {
        sforce.timesis(maxforce/fl);
        // possibly raise an event here
      }
    }
  }
  d->applyforce(&pd,d,&sforce);
  if (g_is_dyna) {
    DL_vector counterforce;
    sforce.neg(&counterforce);
    ((DL_dyna*)g)->applyforce(&pg,g,&counterforce);
  }
}

void DL_spring::springforce(DL_vector *f){
  f->assign(&sforce);
}

void DL_spring::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=1;
  if (g_is_dyna) tf=2; else tf=1;
}

void DL_spring::get_force_info(int i, DL_actuator_type& at,
			      DL_dyna*& _g, DL_point *_p, DL_vector *_v){
  if (i==0) {
    at=force;
    _g=d;
    _p->assign(&pd);
    springforce(_v);
    return;
  }
  if ((i==1) && (g_is_dyna)) {
    at=force;
    _g=(DL_dyna*)g;
    _p->assign(&pg);
    springforce(_v);
    _v->neg(_v);
    return;
  }
  at=none;
};
