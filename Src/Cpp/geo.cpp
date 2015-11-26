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
// filename     : geo.cpp
// description  : non-inline methods of class DL_geo
// author       : Bart Barenbrug     Januari '97
//
 
#include "geo.h"
#include "dyna_system.h"
 
// ************************** //
// non-inline member fuctions //
// ************************** //
 
void DL_geo::move(DL_point *newpos, DL_matrix *neworient){
  DL_matrix Ad;
  DL_vector vtmp;
  DL_Scalar h=DL_dsystem->get_integrator()->stepsize();
  DL_Scalar oldh=DL_dsystem->get_integrator()->old_stepsize();

  if (h==oldh) {
    newpos->minus(&(mstate.z),&(mstate.v));
    mstate.v.timesis(0.5/h);
  }
  else {
    newpos->minus(&(nextmstate.z),&(mstate.v));
    mstate.v.timesis(oldh/((h+oldh)*h));
    nextmstate.z.minus(&(mstate.z),&vtmp);
    vtmp.timesis(h/((h+oldh)*oldh));
    mstate.v.plusis(&vtmp);
  }
  mstate.z.assign(&nextmstate.z);
  nextmstate.z.assign(newpos);
  nextmstate.z.minus(&(mstate.z),&vtmp);
  vtmp.timesis(2.0/h);
  vtmp.minus(&(mstate.v),&(nextmstate.v));

  if (h==oldh) {
    neworient->minus(&(mstate.A),&Ad);
    mstate.A.assign(&(nextmstate.A));
    mstate.w.init(Ad.c0.z*mstate.A.c0.y+Ad.c1.z*mstate.A.c1.y+Ad.c2.z*mstate.A.c2.y,
	  	  Ad.c0.x*mstate.A.c0.z+Ad.c1.x*mstate.A.c1.z+Ad.c2.x*mstate.A.c2.z,
		  Ad.c0.y*mstate.A.c0.x+Ad.c1.y*mstate.A.c1.x+Ad.c2.y*mstate.A.c2.x);
    mstate.w.timesis(0.5/h);
  }
  else {
    nextmstate.A.minus(&(mstate.A),&(mstate.A));
    mstate.A.timesis(h/((h+oldh)*oldh));
    neworient->minus(&(nextmstate.A),&Ad);
    Ad.timesis(oldh/((h+oldh)*h));
    Ad.plusis(&(mstate.A));
    mstate.A.assign(&(nextmstate.A));
    mstate.w.init(Ad.c0.z*mstate.A.c0.y+Ad.c1.z*mstate.A.c1.y+Ad.c2.z*mstate.A.c2.y,
	  	  Ad.c0.x*mstate.A.c0.z+Ad.c1.x*mstate.A.c1.z+Ad.c2.x*mstate.A.c2.z,
		  Ad.c0.y*mstate.A.c0.x+Ad.c1.y*mstate.A.c1.x+Ad.c2.y*mstate.A.c2.x);
  }
  nextmstate.A.assign(neworient);
  nextmstate.A.minus(&(mstate.A),&Ad);
  vtmp.init(Ad.c0.z*neworient->c0.y+Ad.c1.z*neworient->c1.y+Ad.c2.z*neworient->c2.y,
	    Ad.c0.x*neworient->c0.z+Ad.c1.x*neworient->c1.z+Ad.c2.x*neworient->c2.z,
            Ad.c0.y*neworient->c0.x+Ad.c1.y*neworient->c1.x+Ad.c2.y*neworient->c2.x);
  vtmp.timesis(2.0/h);
  vtmp.minus(&(mstate.w),&(nextmstate.w));
}
