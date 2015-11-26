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
// filename: force_drawable.h
// description: class describing entities of which forces can be draw: this
//              provides the interface to obtain the info on those forces.
// author: Bart Barenbrug   Januari '98
//

#ifndef DL_FORCE_DRAWABLEH
#define DL_FORCE_DRAWABLEH

#include "boolean.h"
#include "pointvector.h"

class DL_dyna;

enum DL_actuator_type {none, force, torque, impulse};

// *********************** //
// class DL_force_drawable //
// *********************** //

class DL_force_drawable {
  protected:
    boolean showing;
  public:
    // client data for the force_drawer to adminstrate which object is
    // managing this force_drawable
    void *fd_elem;

    // for (un)registering with the force_drawer:
    virtual void show_forces(void);
    virtual void hide_forces(void);

    // for use by the force_drawer:
    virtual void get_fd_info(int&,int&);
    virtual void get_force_info(int, DL_actuator_type&,
			        DL_dyna*&, DL_point*, DL_vector*);

       DL_force_drawable();  // constructor
       ~DL_force_drawable(); // destructor

};

#endif
