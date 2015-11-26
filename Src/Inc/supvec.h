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
// filename: supvec.h
// description: supvec vectors combine three vectors, one vector4 and one
//              matrix to the motion state of a geo.
//              The operations are dived in two groups: the ones that act
//              on the whole vector (init and assign), and the ones used
//              by the motion integrators of dyna (which only act on the
//              orrientation)
//              The class is internal to the DL module
// author: Bart Barenbrug   May '96


#ifndef DL_SUPVECH
#define DL_SUPVECH

#include <iostream.h>
#include <math.h>
#include "pointvector.h"
#include "matrix.h"
#include "vector4.h"
class DL_dyna;

// *************** //
// class DL_supvec //
// *************** //

class DL_supvec {
  public:   
    DL_point   z;
    DL_vector  v;
    DL_vector4 q;
    DL_vector  w;
    DL_matrix  A;
    
    void       init(DL_point*,DL_vector*,DL_vector4*,DL_vector*, DL_matrix*);
    void       assign(DL_supvec*);

    void       plusis(DL_supvec*);
    void       timesis(DL_Scalar);
    void       times(DL_Scalar,DL_supvec*);       
    void       plus(DL_supvec*,DL_supvec*);

    void       q2A(DL_dyna*);
    void       A2q();
    
    DL_supvec();                                        // constructor
    DL_supvec(DL_point*,DL_vector*,DL_vector4*,DL_vector*,DL_matrix*);
                                                        // constructor
    ~DL_supvec(){};                                     // destructor
};

// handy for debugging:
#define PRINTSUPVEC(S) \
  Msg("z=(%f,%f,%f)\nv=(%f,%f,%f)\nw=(%f,%f,%f)\nq=(%f,%f,%f,%f)\n  (%f %f %f)\nA=(%f %f %f)\n  (%f %f %f)\n", \
      S->z.x, S->z.y, S->z.z, \
      S->v.x, S->v.y, S->v.z, \
      S->w.x, S->w.y, S->w.z, \
      S->q.c[0], S->q.c[1], S->q.c[2], S->q.c[3], \
      S->A.c0->x, S->A.c1->x, S->A.c2->x, \
      S->A.c0->y, S->A.c1->y, S->A.c2->y, \
      S->A.c0->z, S->A.c1->z, S->A.c2->z );
#endif
