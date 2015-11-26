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
// filename: rungekutta.h
// description: 4-point Runge Kutta motion integrator
// author: Bart Barenbrug   March '96
//

#ifndef DL_RUNGEKUTTA4H
#define DL_RUNGEKUTTA4H

#include "m_integrator.h"

// ******************** //
// class DL_rungekutta4 //
// ******************** //

class DL_rungekutta4 : public DL_m_integrator {
  public:
    /// for external (to DL) use:
    DL_Scalar ast(){return 0.5;};

    /// for internal (DL) use only:
    void integrate(DL_supvec*, DL_dyna*, DL_supvec*);
                                 // do one integration step
};

#endif
