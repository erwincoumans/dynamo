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
// filename: csplinesegment.h
// description: cspline segment curve (c-spline with exactly four control
//              points). Used in the construction of c-spline curves with
//              any number of points greater than two.
// author: Bart Barenbrug   April '98
//

#ifndef DL_CSPLINESEGMENTH
#define DL_CSPLINESEGMENTH

#include "curve.h"

// *********************** //
// class DL_csplinesegment //
// *********************** //

class DL_csplinesegment : public DL_curve {
protected:
  DL_point p[4];    // control points
  DL_point a,b,c,d; // coefficients of the powers of the curve parameter.
  DL_Scalar _t,_dt,dtinv;
       // actual curve parameter for this curve is calculated from s as
       // (s-t)/dt (=(s-t)*dtinv).
       // for use by cspline
  void recalc_abcd();
public:
  void assign(DL_csplinesegment*);
       // assignment
  virtual boolean pos(DL_Scalar,DL_point*);
       // returns curve(s) in the point* (in local coordinates of g)
       // and whether s is within bounds as return value
  virtual boolean deriv(DL_Scalar,DL_vector*);
       // returns curve'(s) in the vector* (in local coordinates of g)
       // and whether s in within bounds as return value
          boolean deriv2(DL_Scalar,DL_vector*);
       // returns curve''(s) in the vector* (in local coordinates of g)
       // and whether s in within bounds as return value
  virtual boolean indomain(DL_Scalar);
       // returns whether s is within bounds
  virtual DL_Scalar closeto(DL_point*);
       // returns a curveparameter s with p-curve(s) minimal
       // p given in world coordinates
  void init(DL_geo*,DL_point*, DL_point*,DL_point*,DL_point*);
       // initialise with geo and the four control points (each
       // given in the local coordinates of g)
  void set_interval(DL_Scalar,DL_Scalar);
       // set _t, _dt, minparam and maxparam
  void update_control_point(int,DL_point*);
       // change the coordinates of the i-th control point
  DL_Scalar t() { return _t; }
  DL_Scalar dt() { return _dt; }
  
             DL_csplinesegment();     // constructor
	     ~DL_csplinesegment();    // destructor
};

#endif
