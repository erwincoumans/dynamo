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
// description: class using which a user can draw user forces: this class
//              provides the storage of the required data.
//              Much like gdp's fd_once class, but those forces are only
//              drawn for one frame, with all subsequent creation/destruction
//              of topologies when forces have to be drawn longer.
// author: Bart Barenbrug   Januari '98
//

#ifndef DL_USR_FORCE_DRAWABLEH
#define DL_USR_FORCE_DRAWABLEH

#include "force_drawable.h"

// *************************** //
// class DL_usr_force_drawable //
// *************************** //

class DL_usr_force_drawable : public DL_force_drawable {
  public:
    boolean reaction;
    DL_actuator_type at;
    DL_dyna *d;
    DL_point p;
    DL_vector f;
  
    // for use by the force_drawer:
    virtual void get_fd_info(int&,int&);
    virtual void get_force_info(int, DL_actuator_type&,
			        DL_dyna*&, DL_point*, DL_vector*);

       DL_usr_force_drawable();  // constructor
       DL_usr_force_drawable(boolean react, DL_actuator_type _at,
	    DL_dyna *_d, DL_point *_p, DL_vector *_f);  // constructor

       ~DL_usr_force_drawable(){}; // destructor

};

#endif
