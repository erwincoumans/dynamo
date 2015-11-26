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
// filename     : actuator.cpp
// description	: methods of class DL_actuator
// author	: Bart Barenbrug     March 1998
//

#include "actuator.h"

// *************** //
// member fuctions //
// *************** //

DL_actuator::DL_actuator() : DL_force_drawable() {
  max_actuator=0;
}

DL_actuator::~DL_actuator() {
}

DL_Scalar DL_actuator::limit(DL_Scalar a){
  if (max_actuator>0) {
    if (fabs(a)>max_actuator)
      return (a>0 ? max_actuator : -max_actuator);
  }
  return a;
}

void DL_actuator::apply(DL_Scalar a){
  // nothing to do for the base class: to be implemented by the specialisation
}
