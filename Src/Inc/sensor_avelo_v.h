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
// filename: sensor_avelo_v.h
// description: sensor that measures the angular velocity between two vectors
//              along a given direction vector
// author: Bart Barenbrug   March '98

#ifndef DL_SENSOR_AVELO_VH
#define DL_SENSOR_AVELO_VH

#include "sensor.h"

// *********************** //
// class DL_sensor_avelo_v //
// *********************** //

class DL_sensor_avelo_v: public DL_sensor {
protected:
  DL_geo  *g;
  DL_dyna *d;
  DL_vector vd,rd,vg;
  // the sensor measures the angular velocity between vg in g and vd in d
  // in the direction of rd (vd, vg,rd in local coordinates);
  DL_vector vdw,rdw,vgw;
  // world coordinates of pd,pg and rd;

  boolean g_is_dyna;

  DL_Scalar anglem; // the current angle (mod 360 degrees -180<=anglem<180)
  DL_Scalar d_ang; // angular velocity
  int last_fn; // last frame number
public:
  void init(DL_dyna*, DL_vector*, DL_vector*,
            DL_geo*, DL_vector*);
  void initw(DL_dyna*, DL_vector*, DL_vector*,
             DL_geo*, DL_vector*);
  
  virtual DL_Scalar sense(void);

             DL_sensor_avelo_v();     // constructor
	     ~DL_sensor_avelo_v();    // destructor

};

#endif
