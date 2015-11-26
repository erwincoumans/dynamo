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
// filename: sensor_velo_v.h
// description: sensor that measures the relative velocity between two points
//              along a given direction vector
// author: Bart Barenbrug   March '98

#ifndef DL_SENSOR_VELO_VH
#define DL_SENSOR_VELO_VH

#include "sensor.h"

// ********************** //
// class DL_sensor_velo_v //
// ********************** //

class DL_sensor_velo_v: public DL_sensor {
protected:
  DL_geo  *g;
  DL_dyna *d;
  DL_point pd,pg;
  DL_vector rd;
  // the sensor measures the relative velocity of pg in g and pd in d
  // in the direction of rd (pd, pg,rd in local coordinates);
  DL_vector vdw,vgw;
  // velocity of pd and pw;
  DL_vector rdw;
  // world coordinates of rd;

  boolean g_is_dyna;
public:
  void init(DL_dyna*, DL_point*, DL_vector*,
            DL_geo*, DL_point*);
  
  virtual DL_Scalar sense(void);

             DL_sensor_velo_v();     // constructor
	     ~DL_sensor_velo_v();    // destructor

};

#endif
