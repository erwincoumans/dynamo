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
// filename: connector.h
// description: connector constraint class
// author: Bart Barenbrug   August '96
//

#ifndef DL_CONNECTORH
#define DL_CONNECTORH

#include "constraint_manager.h"
#include "ptp.h"
#include "orientation.h"

// ****************** //
// class DL_connector //
// ****************** //

class DL_connector : public DL_constraint {
public:
  DL_ptp *myptp;
  DL_orientation *myorient;
  void init(DL_dyna*, DL_point*, DL_point*, DL_point*,
            DL_geo* , DL_point*, DL_point*, DL_point*);
	     // initialise the constraint

  // methods for analytical dC/dR determination:
   // the next few methods don't need to do anything: just override
   // them anyway to prevent the errormessage from the constraint-
   // implementations of them...
  virtual boolean dCdRsub(DL_constraint*, DL_largematrix*){return FALSE;}
  virtual boolean dCdFq(DL_dyna*,DL_point*,DL_largematrix*){return FALSE;}
  virtual boolean dCdF(DL_dyna*,DL_largematrix*){return FALSE;}
  virtual boolean dCdM(DL_dyna*,DL_largematrix*){return FALSE;}
  virtual boolean dCdI(DL_dyna*,DL_point*,DL_largematrix*){return FALSE;}
  virtual void apply_restrictions(DL_largevector*){};

  virtual boolean check_restrictions();
                     // check the reaction forces and torques
	             // and maybe deactivate
  virtual void test_restriction_changes(DL_largevector*);
                     // announce to the constraint which restriction change
		     // is about to be applied, so the constraint can
		     // decide to deactivate itself
  // for NOT (un)registering with the force_drawer:
  virtual void show_forces(void){};
  virtual void hide_forces(void){};

             DL_connector();                    // constructor
	     ~DL_connector();                   // destructor
};

#endif
