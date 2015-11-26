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
// filename     : force_drawable.cpp
// description	: non-inline methods of class DL_force_drawable
// author	: Bart Barenbrug     Januari '98
//

#include "force_drawable.h"
#include "force_drawer.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_force_drawable::DL_force_drawable(){
  fd_elem=NULL;
  // don't show forces by default:
  showing=FALSE;
}

DL_force_drawable::~DL_force_drawable() {
  hide_forces();
}

void DL_force_drawable::show_forces(void) {
  if (showing || (!DL_dsystem->get_companion())) return;
  else {
    showing=TRUE;
    DL_dsystem->get_companion()->register_fd(this);
  }
}

void DL_force_drawable::hide_forces(void) {
  if (showing) {
    showing=FALSE;
    DL_dsystem->get_companion()->remove_fd(this);
  }
}

void DL_force_drawable::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  nrf=tf=0;
}

void DL_force_drawable::get_force_info(int i, DL_actuator_type& at,
			      DL_dyna*& g, DL_point *p, DL_vector *v){

  at=none;
};
