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

#include "usr_force_drawable.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_usr_force_drawable::DL_usr_force_drawable():
  DL_force_drawable() {
  reaction=FALSE;
  at=none;
  d=NULL;
  p.init(0,0,0);
  f.init(0,0,0);
}

DL_usr_force_drawable::DL_usr_force_drawable(boolean react,
	DL_actuator_type _at, DL_dyna *_d, DL_point *_p, DL_vector *_f)
      :DL_force_drawable() {
  reaction=react;
  at=_at;
  d=_d;
  p.assign(_p);
  f.assign(_f);
  show_forces();
}

void DL_usr_force_drawable::get_fd_info(int& nrf, int& tf) {
  // forces from 0 to nrf-1;
  // reactionforces from nrf to tf-1;
  if (at==none) nrf=tf=0;
  else {
    nrf=(reaction?0:1);
    tf=1;
  }
}

void DL_usr_force_drawable::get_force_info(int i, DL_actuator_type& _at,
			      DL_dyna*& _d, DL_point* _p, DL_vector* _f){
  if ((at==none) || (i!=0)) _at=none;
  else {
    _at=at;
    _d=d;
    _p->assign(&p);
    _f->assign(&f);
  }
};
