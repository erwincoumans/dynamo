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
// filename     : controller.cpp
// description	: non-inline methods of class DL_controller
// author	: Bart Barenbrug     Februari '97
//

#include "controller.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_controller::DL_controller(): DL_ListElem(), DL_force_drawable(){
  active=FALSE;
}

DL_controller::~DL_controller() {
  deactivate();
}

void DL_controller::activate(void) {
  if (active) return;
  else {
    DL_dsystem->add_controller(this);
    active=TRUE;
  }
}

void DL_controller::deactivate(void) {
  if (active) {
    DL_dsystem->rem_controller(this);
    active=FALSE;
  }
}

