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
// filename     : m_integrator.cpp
// description	: non-inline methods of class m_integrator
// author	: Bart Barenbrug     March '96
//

#include "m_integrator.h"
#include "dyna.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_Scalar DL_m_integrator::stepsize() {
  return h;
}

DL_Scalar DL_m_integrator::halfstepsize() {
  return halfh;
}

void DL_m_integrator::set_stepsize(DL_Scalar newh) {
  h=newh; halfh=0.5*newh;
}
