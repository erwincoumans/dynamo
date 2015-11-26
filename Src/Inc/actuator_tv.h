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
// filename: actuator_tv.h
// description: actuator which applies a torque along a given vector
// author: Bart Barenbrug   March '98

#ifndef DL_ACTUATOR_TVH
#define DL_ACTUATOR_TVH

#include "actuator.h"

// ******************** //
// class DL_actuator_tv //
// ******************** //

class DL_actuator_tv: public DL_actuator {
protected:
  DL_geo  *g;
  DL_dyna *d;
  DL_vector rd;
  // the actuator applies a torque in the direction of rd to pg in g
  // and pd in d (pd, pg,rd in local coordinates);
  DL_vector t;
  // calculated torque

  boolean g_is_dyna;
  
public:
  virtual void apply(DL_Scalar);

  void init(DL_dyna*, DL_vector*, DL_geo*);
  void initw(DL_dyna*, DL_vector*, DL_geo*);
  
            DL_actuator_tv();    // constructor
	    ~DL_actuator_tv();   // destructor

  virtual void get_fd_info(int&,int&);
  virtual void get_force_info(int, DL_actuator_type&,
		              DL_dyna*&, DL_point*, DL_vector*);

  
};

#endif
