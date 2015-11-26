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

// filename    : matrix.h

#ifndef DL_MATRIXH
#define DL_MATRIXH

#include <iostream.h>
#include "pointvector.h"

class DL_matrix { 
  public:  
    DL_vector c0;         // colomn 0
    DL_vector c1;         // colomn 1
    DL_vector c2;         // colomn 2

    void  makeone();
    void  makezero();
    void  normalize();
    
    void  assign(DL_matrix*);
    void  assign(DL_vector *v0,DL_vector *v1,DL_vector *v2);
    void  assign(DL_Scalar c0x, DL_Scalar c1x, DL_Scalar c2x,
                 DL_Scalar c0y, DL_Scalar c1y, DL_Scalar c2y,
		 DL_Scalar c0z, DL_Scalar c1z, DL_Scalar c2z);
    
    DL_Scalar get(int,int);
    void  set(int,int,DL_Scalar);
    void  plus(DL_matrix*,DL_matrix*);
    void  minus(DL_matrix*,DL_matrix*);

    void  plusis(DL_matrix*);
    void  minusis(DL_matrix*);
    void  timesis(DL_Scalar);
    
    void  times(DL_matrix*,DL_matrix*);
    void  times(DL_Scalar,DL_matrix*);
    void  times(DL_vector*,DL_vector*);
    void  times(DL_point*,DL_point*);
    void  invert(DL_matrix*);
    void  transpose(DL_matrix*);
    void  timestranspose(DL_matrix*, DL_matrix*);
    void  transposetimes(DL_vector*, DL_vector*);
    void  jacobi(DL_matrix*,DL_vector*);
    void  diag_transpose_vec(DL_vector*,DL_vector*,DL_vector*);
    void  negcrossdiagcross(DL_point*,DL_vector*,DL_point*);
    void  diagcrosstranspose(DL_vector*,DL_point*,DL_matrix*);

    void  tensor(DL_vector*,DL_vector*);

    DL_matrix(){};
    DL_matrix(DL_vector* d0, DL_vector* d1, DL_vector* d2);
    ~DL_matrix(){};

};

#endif
