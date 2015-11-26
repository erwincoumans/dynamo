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
// filename: vector4.h
// description: quaternions (required by forward dynamics for
//              orientation storage)
//              This class is internal to the DL module
// author: Bart Barenbrug   March '96
//

#ifndef DL_VECTOR4H
#define DL_VECTOR4H

#include <math.h>
#include "pointvector.h"
#include "matrix.h"

// **************** //
// class DL_vector4 //
// **************** //

class DL_vector4 {
  public:   
    DL_Scalar  c[4];
    void       init(DL_Scalar,DL_Scalar,DL_Scalar,DL_Scalar);
    void       assign(DL_vector4*);
    DL_Scalar  norm();
    void       normalize();
    DL_Scalar  inprod(DL_vector4*);
    void       plusis(DL_vector4*);
    void       timesis(DL_Scalar);
    void       timesis(DL_vector4*);
    boolean    equal(DL_vector4*);

    DL_Scalar  get(int);
    void       set(int,DL_Scalar);

    void       neg(DL_vector4*);
    void       times(DL_Scalar,DL_vector4*);       
    void       plus(DL_vector4*,DL_vector4*);      
    void       minus(DL_vector4*,DL_vector4*);     
    void       times(DL_vector*,DL_vector4*);

    // methods that assume this vector represents an orientation
    // in Euler coordinates:
    void       to_matrix(DL_matrix*);
    void       from_matrix(DL_matrix*);

               DL_vector4(){};                       // constructor        
               DL_vector4(DL_Scalar,DL_Scalar,DL_Scalar,DL_Scalar);
                                                     // constructor
	       ~DL_vector4(){};                      // destructor
};

#endif
