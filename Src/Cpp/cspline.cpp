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
// filename     : cspline.cpp
// description	: non-inline methods of class DL_cspline
// author	: Bart Barenbrug     April '98
//

#include "cspline.h"
#include "dyna_system.h"

// ************************** //
// non-inline member fuctions //
// ************************** //

DL_cspline::DL_cspline():DL_curve() {
  cyclic=FALSE;
  segment=NULL;
}

DL_cspline::~DL_cspline(void) {
  // delete all segments:
  if (segment) delete[] segment;
}

void DL_cspline::init(DL_geo *geo, DL_List *points, boolean _cyclic) {
  DL_point *p0,*p1,*p2,*p3;
  int i,n=points->length();
  
  if (n<3) {
    DL_dsystem->get_companion()->Msg("Error: cspline::init: a cspline requires at least 3 control points\n cspline curve not initialised\n");
    return;
  }
  DL_curve::init(geo);
  cyclic=_cyclic;

  minparam=0;
  maxparam=(float)( cyclic ? n : n-1 );
  
  // initialise all segments;
  segment=new DL_csplinesegment[(int)maxparam];
  
  if (cyclic) {
    p0=(DL_point*)(points->getfirst());
    p1=(DL_point*)(points->getnext(p0));
    p2=(DL_point*)(points->getnext(p1));
    p3=(DL_point*)(points->getnext(p2));
    for (i=0; i<n; i++) {
      segment[i].init(g,p0,p1,p2,p3);
      segment[i].set_interval(i,1); // set t and dt
      p0=p1; p1=p2; p2=p3; p3=(DL_point*)(points->getnext(p2));
      if (p3==NULL) p3=(DL_point*)(points->getfirst());
    }
  }
  else {
    p1=p2=(DL_point*)(points->getfirst());
    p3=(DL_point*)(points->getnext(p2));
    for (i=0; i<n-1; i++) {
      p0=p1; p1=p2; p2=p3; p3=(DL_point*)(points->getnext(p2));
      if (p3==NULL) p3=p2;
      segment[i].init(g,p0,p1,p2,p3);
      segment[i].set_interval(i,1); // set t and dt
    }
  }
}

void DL_cspline::assign(DL_cspline* bs){
  if (bs) {
    // now we have to make a deep copy of the list:
    // first delete the old list and make an empty new one:
    int i;
    delete[] segment;
    segment=new DL_csplinesegment[(int)bs->maxparam];

    // then copy:
    for (i=0; i<(int)bs->maxparam; i++) {
      segment[i].assign(&(bs->segment[i]));
    }

    DL_curve::assign(bs);
    cyclic=bs->cyclic;
  }
}

boolean DL_cspline::pos(DL_Scalar s, DL_point *p) {
// returns curve(s) in the point* (in local coordinates of g) and
// whether s is within bounds as return value
  int i; // index specifying correct segment
  if (cyclic) {
    int sn=(int)floor(s);
    i=sn%(int)maxparam;
    if (i<0) i+=(int)maxparam;
    s-=(sn-i);
  }
  else {
    if (s<0) {
      segment[0].pos(0,p);
      return FALSE;
    }
    if (s>maxparam) {
      segment[(int)maxparam-1].pos(maxparam,p);
      return FALSE;
    }
    i=(int)floor(s);
    if (i==(int)maxparam) i-=1;
  }
  // ok: we're in bounds
  segment[i].pos(s,p);
  return TRUE;
}

boolean DL_cspline::deriv(DL_Scalar s, DL_vector *v){
// returns curve'(s) in the vector* (in local coordinates of g)
// and whether s in within bounds as return value
  int i; // index specifying correct segment
  if (cyclic) {
    int sn=(int)floor(s);
    i=sn%(int)maxparam;
    if (i<0) i+=(int)maxparam;
    s-=(sn-i);
  }
  else {
    if (s<0) {
      segment[0].deriv(0,v);
      return FALSE;
    }
    if (s>maxparam) {
      segment[(int)maxparam].deriv(maxparam,v);
      return FALSE;
    }
    i=(int)floor(s) % (int)maxparam;
  }
  // ok: we're in bounds
  segment[i].deriv(s,v);
  return TRUE;
}

boolean DL_cspline::deriv2(DL_Scalar s, DL_vector *v){
// returns curve''(s) in the vector* (in local coordinates of g)
// and whether s in within bounds as return value
  int i; // index specifying correct segment
  if (cyclic) {
    int sn=(int)floor(s);
    i=sn%(int)maxparam;
    if (i<0) i+=(int)maxparam;
    s-=(sn-i);
  }
  else {
    if (s<0) {
      segment[0].deriv2(0,v);
      return FALSE;
    }
    if (s>maxparam) {
      segment[(int)maxparam].deriv2(maxparam,v);
      return FALSE;
    }
    i=(int)floor(s) % (int)maxparam;
  }
  // ok: we're in bounds
  segment[i].deriv2(s,v);
  return TRUE;
}

boolean DL_cspline::indomain(DL_Scalar s) {
// returns whether s is within bounds
  return (cyclic || ((0<=s) && (s<=maxparam)));
}

DL_Scalar DL_cspline::closeto(DL_point *p) {
// returns a curveparameter s with p-curve(s) minimal
// p given in world coordinates
#define STEPSIZE 0.04
  // b-splines are rather smooth, so we just sample the curve.
  DL_Scalar s,s_min,mind=100000,dist;
  DL_point pl,pc;
  DL_vector pdiff,d;
  int i;

  // first convert p to local coordinates
  if (g) g->to_local(p,NULL,&pl);
  else pl.assign(p);

  s=0;
  while (s<=maxparam) {
    pos(s,&pc);
    pc.minus(&pl,&pdiff);
    if ((dist=pdiff.norm())<mind) {
      s_min=s; mind=dist;
    }
    s+=STEPSIZE;  // assume equal dt's for all segments
  }
#undef STEPSIZE

  // now we're pretty close. Do a local search for refinement:
  // maximum 10 steps, or when there's not enough movement anymore. 
  s=s_min;
  i=0;
  do {
    pos(s,&pc);
    deriv(s,&d);
    pl.minus(&pc,&pdiff);
    s+=(dist=pdiff.inprod(&d)/d.inprod(&d));
    i++;
  } while ((i<10) && (fabs(dist)<0.0001));
  
  return s;
}

void DL_cspline::update_control_point(int i, DL_point *pn){
  int j,segm,M=(int)maxparam;
  if ((i<0) || (i> (cyclic ? M-1 : M))) {
    DL_dsystem->get_companion()->Msg("Warning: DL_cspline::update_control_point(int,DL_point): index out of range\n");
    return;
  }
  if (cyclic) {
    for (j=0;j<4;j++) {
      segm=i-j; if (segm<0) segm+=M;
      segment[segm].update_control_point(j,pn);
    }
  }
  else {
    if (i==0) {
      segment[0].update_control_point(1,pn);
      segment[0].update_control_point(0,pn);
      segment[1].update_control_point(0,pn);
    }
    if (i==1) {
      segment[0].update_control_point(2,pn);
      segment[1].update_control_point(1,pn);
      if (maxparam>2) segment[2].update_control_point(0,pn);
    } else
    if (i+1==M) {
      segment[M-3].update_control_point(3,pn);
      segment[M-2].update_control_point(2,pn);
      segment[M-1].update_control_point(1,pn);           
    } else
    if (i==M) {
      segment[M-2].update_control_point(3,pn);
      segment[M-1].update_control_point(3,pn);
      segment[M-1].update_control_point(2,pn);
    }
    else for (j=0;j<4;j++) segment[i+1-j].update_control_point(j,pn);
  }
}
