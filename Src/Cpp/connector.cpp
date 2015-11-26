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
// filename     : connector.cpp
// description	: non-inline methods of class connector
// author	: Bart Barenbrug     August '96
//

#include "dyna_system.h"
#include "connector.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_connector::DL_connector():DL_constraint() {
  dim=0;
  myorient=new DL_orientation;
  myptp=new DL_ptp;
}

DL_connector::~DL_connector() {
  delete myorient; delete myptp;
}

void DL_connector::init(
    DL_dyna* _d, DL_point* _pd0, DL_point *_pd1, DL_point *_pd2,
    DL_geo* _g,  DL_point* _pg0, DL_point *_pg1, DL_point *_pg2) {
  if (_d==_g) {
    DL_dsystem->get_companion()->Msg("Error: connector::init: a connector constraint needs points from _different_ objects\n connector constraint not initialised\n");
    return;
  }

  myptp->init(_d,_pd0,_g,_pg0);

  DL_vector v0,v1,v2,w0,w1,w2;

  _pd1->minus(_pd0,&v0);
  _pd2->minus(_pd0,&v1);
  v0.crossprod(&v1,&v2);
  v0.crossprod(&v2,&v1);
  
  _pg1->minus(_pg0,&w0);
  _pg2->minus(_pg0,&w1);
  w0.crossprod(&w1,&w2);
  w0.crossprod(&w2,&w1);

  myorient->init(_d,&v0,&v1,_g,&w1,&w2);
    
  DL_constraint::init();
}

void DL_connector::test_restriction_changes(DL_largevector* lv) {
  if (!myorient->active) {
    myptp->deactivate();
    deactivate();
    // possibly raise an event here
    DL_dsystem->get_companion()->Msg(" connector constraint deactivated\n");
    return;
  }
  if (!myptp->active) {
    myorient->deactivate();
    deactivate();
    // possibly raise an event here
    DL_dsystem->get_companion()->Msg(" connector constraint deactivated\n");
    return;
  }
}

boolean DL_connector::check_restrictions() {
  if (!myorient->active) {
    myptp->deactivate();
    deactivate();
    // possibly raise an event here
    DL_dsystem->get_companion()->Msg(" connector constraint deactivated\n");
    return FALSE;
  }
  if (!myptp->active) {
    myorient->deactivate();
    deactivate();
    // possibly raise an event here
    DL_dsystem->get_companion()->Msg(" connector constraint deactivated\n");
    return FALSE;
  }
  return TRUE;
}
