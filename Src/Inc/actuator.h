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
// filename: actuator.h
// description: base class for making actuators (used by controllers).
// author: Bart Barenbrug   March '98

#ifndef DL_ACTUATORH
#define DL_ACTUATORH

#include "force_drawable.h"

// ***************** //
// class DL_actuator //
// ***************** //

class DL_actuator : public DL_force_drawable {
protected:
  DL_Scalar max_actuator;
  DL_Scalar limit(DL_Scalar); // limit the input against max_actuator
public:
  /// for use by the controllers:
  virtual void apply(DL_Scalar);

  /// for external use:
  void set_max_actuator(DL_Scalar ma){ max_actuator=ma; };
  DL_Scalar get_max_actuator() { return max_actuator; };
  
             DL_actuator();     // constructor
	     ~DL_actuator();    // destructor
};

#endif
