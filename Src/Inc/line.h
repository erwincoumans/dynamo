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
// filename: straightcurve.h
// description: a line curve class (for use with
//              point-to-curve constraint)
// author: Bart Barenbrug   June '96
//

#ifndef DL_LINEH
#define DL_LINEH

#include "curve.h"

// ************* //
// class DL_line //
// ************* //

class DL_line : public DL_curve {
protected:
  DL_point pnt;  // in local coordinates of geo
  DL_vector dir; // in local coordinates of geo. This is line x: pnt+x*dir
                 // dir is normalised
public:
  void assign(DL_line*);
       // assignment
  virtual boolean pos(DL_Scalar,DL_point*);
       // returns curve(s) in the point* (in local coordinates of g)
       // and whether s is within bounds as return value
  virtual boolean deriv(DL_Scalar,DL_vector*);
       // returns curve'(s) in the vector* (in local coordinates of g)
       // and whether s in within bounds as return value
  virtual boolean indomain(DL_Scalar);
       // returns whether s is within bounds
  virtual DL_Scalar closeto(DL_point*);
       // returns a curveparameter s with p-curve(s) minimal
  void init(DL_geo*, DL_point*, DL_vector*); // initialise with geo, p and d
  void set_minparam(DL_Scalar);
  void set_maxparam(DL_Scalar);
  
             DL_line();     // constructor
	     ~DL_line();    // destructor
};

#endif
