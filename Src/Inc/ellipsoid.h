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
// Filename: ellipsoid.h
// description: ellipsoid surface class (for use with
//              point-to-surface constraint)
// author: Bart Barenbrug   July '96
//

#ifndef DL_ELIPSOIDH
#define DL_ELIPSOIDH

#include "surface.h"

// ***************** //
// class DL_ellipsoid //
// ***************** //

class DL_ellipsoid : public DL_surface {
protected:
  // this is the ellipsoid with center c and major axis x, y and z
  DL_point c;   // in local coordinates of geo
  DL_vector x, y, z; // in local coordinates of geo
public:
  void assign(DL_ellipsoid*);
       // assignment
  virtual boolean pos(DL_Scalar,DL_Scalar,DL_point*);
       // returns surface(s,t) (in local coordinates of g) in the point* and
       // whether s,t is within bounds as return value
  virtual boolean deriv0(DL_Scalar,DL_Scalar,DL_vector*);
       // returns dsurface(s,t)/ds (in local coordinates of g) in the vector*
       // and whether s,t is within bounds as return value
  virtual boolean deriv1(DL_Scalar,DL_Scalar,DL_vector*);
       // returns dsurface(s,t)/dt (in local coordinates of g) in the vector*
       // and whether s,t is within bounds as return value
  virtual boolean indomain(DL_Scalar,DL_Scalar);
       // returns whether s,t is within bounds
  virtual boolean closeto(DL_point*,DL_Scalar&,DL_Scalar&);
       // calculates curveparameters s,t with p-surface(s,t) minimal
       // returns whether those s,t are within bounds
  void init(DL_geo*,DL_point*,DL_vector*,DL_vector*,DL_vector*);
       // initialise with geo, c,x,y,z

             DL_ellipsoid();     // constructor
	     ~DL_ellipsoid();    // destructor
};

#endif
