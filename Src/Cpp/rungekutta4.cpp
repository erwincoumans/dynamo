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
// filename     : rungekutta4.cpp
// description	: non-inline methods of class rungekutta
// author	: Bart Barenbrug     March '96
//

#include "rungekutta4.h"
#include "dyna.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

void DL_rungekutta4::integrate(DL_supvec *y, DL_dyna *d, DL_supvec *ny) {
  DL_supvec k1, k2, k3, k4;

  d->ode(y,&k1,0);  // dfh==0, so use the constant instead)
  k1.times(halfh,ny);
  ny->plusis(y);
  ny->q2A(d);

  d->ode(ny,&k2,0);
  k2.times(halfh,ny);
  ny->plusis(y);
  ny->q2A(d);
  
  d->ode(ny,&k3,0);
  k3.times(h,ny);
  ny->plusis(y);
  ny->q2A(d);

  d->ode(ny,&k4,0);

  k2.timesis(2);
  k3.timesis(2);

  k1.plus(&k2,ny);
  ny->plusis(&k3);
  ny->plusis(&k4);
  ny->timesis((1.0/6.0)*h);
  ny->plusis(y);

  ny->q2A(d);
}
