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
// filename: rope.h
// description: rope class: controller to switch the bar-constraint
//              on and off if it it used as a rope. Not to be used
//              on its own.
// author: Bart Barenbrug   March '97
//

#ifndef DL_ROPEH
#define DL_ROPEH

#include "controller.h"
#include "dyna_system.h"

// ************* //
// class DL_rope //
// ************* //

// forward declaration:
class DL_bar;

class DL_rope : public DL_controller {
protected:
  DL_bar *b;
  DL_point *pd,*pg;
  DL_Scalar lsqr;
public:
  /// for use by DL_bar only:

       void init(DL_bar*, DL_point*, DL_point*, DL_Scalar);
       void set_lsqr(DL_Scalar l){lsqr=l;}
  
       DL_rope();  // constructor
       ~DL_rope(); // destructor

  /// for (DL) internal use only:
  virtual void calculate_and_apply(void);
                           // calculate and apply the controller
			   // forces/torques (called by the dyna system)
};

#endif

