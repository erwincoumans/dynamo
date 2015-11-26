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
// filename: controller.h
// description: class describing controllers in general. To be used
//              to derive ones own controllers from.
// author: Bart Barenbrug   Februari '97
//

#ifndef DL_CONTROLLERH
#define DL_CONTROLLERH

#include "boolean.h"
#include "list.h"
#include "force_drawable.h"

// ******************* //
// class DL_controller //
// ******************* //

class DL_controller : public DL_ListElem, public DL_force_drawable {
protected:
  boolean active;
public:
  /// for external (to DL) use:
  virtual void activate(void);     // (re-activate) the controller
  virtual void deactivate(void);   // deactivate the controller
  boolean is_active(void){return active; };

       DL_controller();                 // constructor
       ~DL_controller(); // destructor

  /// for (DL) internal use only:
  virtual void calculate_and_apply(void)=0;
                           // calculate and apply the controller
			   // forces/torques (called by the dyna system)
};

#endif
