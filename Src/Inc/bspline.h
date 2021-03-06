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
// filename: bspline.h
// description: bspline curve (consisting of b-spline segments)
//              This class is modelled after the splines from Walt
//              with the simplification that no parameter intervals
//              can be specified by the user (all spline segments
//              cover a parameter interval of 1).
// author: Bart Barenbrug   June '96
//

#ifndef DL_BSPLINEH
#define DL_BSPLINEH

#include "bsplinesegment.h"

// **************** //
// class DL_bspline //
// **************** //

class DL_bspline : public DL_curve {
protected:
  boolean cyclic; // is the curve cyclic?
  DL_bsplinesegment *segment;  // the bspline segments
public:
  void assign(DL_bspline*);
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
  void init(DL_geo*,DL_List*,boolean);
       // initialise with geo, points and _cyclic;
  void update_control_point(int,DL_point*);
       // change the coordinates of the i-th control point

             DL_bspline();     // constructor
	     ~DL_bspline();    // destructor
};

#endif
