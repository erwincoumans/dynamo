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
// filename     : doubleeuler.cpp
// description	: non-inline methods of class DL_double_euler
// author	: Bart Barenbrug     May '96
//

#include "doubleeuler.h"
#include "dyna.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

#define dfactor 0.25 // actually: this is half the dfactor since
                     // double_euler takes step of h/2

DL_double_euler::DL_double_euler():DL_m_integrator() {
  dfh=dfactor*h;
}

void DL_double_euler::set_stepsize(DL_Scalar newh) {
  DL_m_integrator::set_stepsize(newh);
  dfh=dfactor*h;
}

void DL_double_euler::integrate(DL_supvec *y, DL_dyna *d, DL_supvec *ny) {
  DL_supvec k1, k2;

  d->ode(y,&k1,dfh);
  k1.times(halfh,ny);
  ny->plusis(y);
  ny->q2A(d);

  d->ode(ny,&k2,dfh);
  k2.timesis(halfh);
  ny->plusis(&k2);
  ny->q2A(d);
}

#undef dfactor
