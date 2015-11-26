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
// filename: force_drawer.h
// description: interface to DL_forcedrawable for a force_drawer
// author: Bart Barenbrug   Januari '98
//

#ifndef DL_FORCE_DRAWERH
#define DL_FORCE_DRAWERH

#include "force_drawable.h"

// ********************* //
// class DL_force_drawer //
// ********************* //

class DL_force_drawer {
public:
  virtual void draw_force(DL_dyna*,DL_point*,DL_vector*){};
  virtual void draw_torque(DL_dyna*,DL_vector*){};
  virtual void draw_impulse(DL_dyna*,DL_point*,DL_vector*){};

  virtual void draw_reaction_force(DL_dyna*,DL_point*,DL_vector*){};
  virtual void draw_reaction_torque(DL_dyna*,DL_vector*){};
  virtual void draw_reaction_impulse(DL_dyna*,DL_point*,DL_vector*){};

  virtual void register_fd(DL_force_drawable*){};
  virtual void remove_fd(DL_force_drawable*){};

       DL_force_drawer(){};  // constructor
       ~DL_force_drawer(){}; // destructor

};

#endif
