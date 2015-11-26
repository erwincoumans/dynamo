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
// filename     : rope.cpp
// description	: non-inline methods of class DL_rope
// author	: Bart Barenbrug     March '97
//

#include "bar.h"
#include "rope.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_rope::DL_rope(): DL_controller(){
  b=NULL; pd=pg=NULL; lsqr=0;
}

DL_rope::~DL_rope(){
  deactivate();
}

void DL_rope::init(DL_bar *_b, DL_point *_pd, DL_point *_pg, DL_Scalar _lsqr){
  b=_b;
  pd=_pd;
  pg=_pg;
  lsqr=_lsqr;
}

void DL_rope::calculate_and_apply(void){
   DL_vector pdiff;
   DL_point pdw,pgw;

   b->get_dyna()->new_toworld(pd,&pdw);
   if (b->get_geo()) b->get_geo()->new_toworld(pg,&pgw);
   else pgw.assign(pg);
   pdw.minus(&pgw,&pdiff);
   if (pdiff.inprod(&pdiff)>=lsqr) {
     b->activate();
     deactivate();
   }
}

