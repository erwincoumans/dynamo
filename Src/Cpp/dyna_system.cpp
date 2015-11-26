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
// filename     : dyna_system.cpp
// description	: non-inline methods of class DL_dyna_system
// author	: Bart Barenbrug     March '96
//

#include "dyna.h"
#include "dyna_system.h"
#include "constraint_manager.h"

// pointer to the one and only dyna_system:
DL_dyna_system* DL_dsystem=NULL;

// ********************** //
// public member fuctions //
// ********************** //

DL_dyna_system::DL_dyna_system(DL_dyna_system_callbacks* comp,
			       DL_m_integrator *mi){
  if (DL_dsystem && DL_dsystem!=this)
    companion->Msg("Warning: there should only be one dyna_system component!!\n");
  if (!DL_dsystem) DL_dsystem=this;
  companion=comp;
  integrator=mi;
  gravity.init(0,0,0);
  frame_nr=0;
  curtime=0;
  srand(12345);
}

void DL_dyna_system::set_gravity(DL_vector* v) {
  gravity.assign(v);
}

void DL_dyna_system::get_gravity(DL_vector *g) {
  g->assign(&gravity);
}

DL_geo* DL_dyna_system::register_geo(void *comp){
  // traverse the geos-list to see if we already have a geo listed
  // with comp as its companion. If so: return a reference to this geo,
  // if not: create a new geo as companion for comp and add it to the
  // list (and return a reference to it).
  if (!comp) return NULL;
  DL_geo *g=(DL_geo*)geos.getfirst();
  while (g) {
    if (g->get_companion()==comp) return g;
    g=(DL_geo*)geos.getnext(g);
  }
  g=new DL_geo(comp);
  companion->get_first_geo_info(g);
  geos.addelem(g);
  return g;
}

void DL_dyna_system::remove_geo(DL_geo *g){
  geos.remelem(g);
}

void DL_dyna_system::register_dyna(DL_dyna *d){
  dynas.addelem(d);
}

void DL_dyna_system::remove_dyna(DL_dyna *d){
  dynas.remelem(d);
}

void DL_dyna_system::add_controller(DL_controller *c) {
  if (show_con_forces) c->show_forces();
  else c->hide_forces();
  controllers.addelem(c);
}

void DL_dyna_system::rem_controller(DL_controller *c) {
  controllers.remelem(c);
}

void DL_dyna_system::update_dyna_companions(void) {
  DL_dyna *d=(DL_dyna*)dynas.getfirst();
  while (d) {
    companion->update_dyna_companion(d);
    d=(DL_dyna*)dynas.getnext(d);
  }
}

void DL_dyna_system::dynamics(void) {
  DL_controller *con=(DL_controller*)controllers.getfirst();
  while(con) {
    con->calculate_and_apply();
    con=(DL_controller*)controllers.getnext(con);
  }
  DL_dyna *d=(DL_dyna*)dynas.getfirst();
  while (d) {
    companion->check_inertiatensor(d);
    d->new_frame();
    d=(DL_dyna*)dynas.getnext(d);
  }
  DL_geo *g=(DL_geo*)geos.getfirst();
  while (g) {
    companion->get_new_geo_info(g);
    g=(DL_geo*)geos.getnext(g);
  }
  
  if (DL_constraints) DL_constraints->satisfy();

  d=(DL_dyna*)dynas.getfirst();
  while (d) {
    d->prepare_for_next_frame();
    d=(DL_dyna*)dynas.getnext(d);
  }
  update_dyna_companions();

  if ((gravity.x!=0.0)||(gravity.y!=0.0)||(gravity.z!=0.0)) {
    DL_vector g;
    d=(DL_dyna*)dynas.getfirst();
    while (d) {
      gravity.times(d->get_mass(),&g);
      d->applycenterforce(&g);
      d=(DL_dyna*)dynas.getnext(d);
    }
  }
  frame_nr++;
  if (integrator) {
    curtime+=integrator->stepsize();
    integrator->shift_stepsize();
  }
}

// C entry-point for the module to call dynamics each frame
extern "C" {
  void DL_dynamics(void) {
    if (DL_dsystem) DL_dsystem->dynamics();
  }
}

void DL_dyna_system::set_integrator(DL_m_integrator *i) {
  if (!i) {
    companion->Msg("Warning: DL_dyna_system::set_integrator(): an integrator is required.\n         Setting integrator to null ignored.\n");
  }
  else {
    if (integrator) { // retain stepsize
      i->set_stepsize(integrator->old_stepsize());
      i->shift_stepsize();
      i->set_stepsize(integrator->stepsize());
    }
    integrator=i;
    // warn all the dyna's that they have to re-integrate
    // using the new integrator
    DL_dyna *d=(DL_dyna*)dynas.getfirst();
    while (d) {
      d->reintegrate();
      d=(DL_dyna*)dynas.getnext(d);
    }
  }
}

DL_m_integrator* DL_dyna_system::get_integrator(void) {
  return(integrator);
}

DL_Scalar DL_dyna_system::kinenergy(void) {
  DL_Scalar tot=0;
  DL_dyna *d=(DL_dyna*)dynas.getfirst();
  while (d) {
    tot+=d->kinenergy();
    d=(DL_dyna*)dynas.getnext(d);
  }
  return tot;
}

DL_Scalar DL_dyna_system::potenergy(void) {
  DL_Scalar tot=0;
  DL_dyna *d=(DL_dyna*)dynas.getfirst();
  while (d) {
    tot+=d->potenergy();
    d=(DL_dyna*)dynas.getnext(d);
  }
  return tot;
}

DL_Scalar DL_dyna_system::totenergy(void) {
  DL_Scalar tot=0;
  DL_dyna *d=(DL_dyna*)dynas.getfirst();
  while (d) {
    tot+=d->totenergy();
    d=(DL_dyna*)dynas.getnext(d);
  }
  return tot;
}

DL_Scalar DL_dyna_system::newkinenergy(void) {
  DL_Scalar tot=0;
  DL_dyna *d=(DL_dyna*)dynas.getfirst();
  while (d) {
    tot+=d->newkinenergy();
    d=(DL_dyna*)dynas.getnext(d);
  }
  return tot;
}

DL_Scalar DL_dyna_system::newpotenergy(void) {
  DL_Scalar tot=0;
  DL_dyna *d=(DL_dyna*)dynas.getfirst();
  while (d) {
    tot+=d->newpotenergy();
    d=(DL_dyna*)dynas.getnext(d);
  }
  return tot;
}

DL_Scalar DL_dyna_system::newtotenergy(void) {
  DL_Scalar tot=0;
  DL_dyna *d=(DL_dyna*)dynas.getfirst();
  while (d) {
    tot+=d->newtotenergy();
    d=(DL_dyna*)dynas.getnext(d);
  }
  return tot;
}

void DL_dyna_system::show_controller_forces() {
  if (show_con_forces) return;
  DL_controller *cc=(DL_controller*)controllers.getfirst();
  while (cc) {
    cc->show_forces();
    cc=(DL_controller*)controllers.getnext(cc);
  }
  show_con_forces=TRUE;
}

void DL_dyna_system::hide_controller_forces() {
  if (!show_con_forces) return;
  DL_controller *cc=(DL_controller*)controllers.getfirst();
  while (cc) {
    cc->hide_forces();
    cc=(DL_controller*)controllers.getnext(cc);
  }
  show_con_forces=FALSE;
}
  
