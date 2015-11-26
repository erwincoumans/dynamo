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
// filename     : sensor_angle_v.cpp
// description	: methods of class DL_sensor_angle_v
// author	: Bart Barenbrug     March 1998
//

#include "dyna.h"
#include "sensor_angle_v.h"

// *************** //
// member fuctions //
// *************** //

DL_sensor_angle_v::DL_sensor_angle_v() : DL_sensor() {
  angle=anglem=0;
}

DL_sensor_angle_v::~DL_sensor_angle_v(void) {
}

void DL_sensor_angle_v::initw(DL_dyna *_d, DL_vector *_vd, DL_vector *_rd,
                              DL_geo  *_g, DL_vector* _vg){
  DL_vector vgl, vdl, rdl;
  _d->to_local(_vd,NULL,&vdl);
  _d->to_local(_rd,NULL,&rdl);
  _g->to_local(_vg,NULL,&vgl);
  init(_d, &vdl, &rdl, _g, &vgl);
}

void DL_sensor_angle_v::init(DL_dyna *_d, DL_vector *_vd, DL_vector *_rd,
                             DL_geo  *_g, DL_vector* _vg){
  if (_g==_d) {
    DL_dsystem->get_companion()->Msg("Error: sensor_angle_v::init: two different geometries are required!\n sensor not initialised\n");
    return;
  }

  d=_d;
  g=_g;
  vd.assign(_vd); vd.normalize();
  vg.assign(_vg); vg.normalize();
  rd.assign(_rd); rd.normalize();
  if ((fabs(vd.norm()-1)>0.001) ||
      (fabs(vg.norm()-1)>0.001) ||
      (fabs(rd.norm()-1)>0.001)) {
    DL_dsystem->get_companion()->Msg("Error: sensor_angle_v::init: vectors should be non-zero\n sensor not initialised\n");
    return;
  }

  // remove any component of vd in the direction of rd:
  DL_vector tmp;
  rd.times(-rd.inprod(&vd),&tmp); vd.plusis(&tmp);
  vd.normalize();
 
  if (g) g_is_dyna=g->is_dyna();
  else g_is_dyna=FALSE;
}

DL_Scalar DL_sensor_angle_v::sense(){
  DL_Scalar ang, d_ang;
  g->to_world(&vg,&vgw);
  d->to_world(&vd,&vdw);
  d->to_world(&rd,&rdw);
  // remove any component of vgw in the direction of rdw:
  DL_vector tmp;
  rdw.times(-rdw.inprod(&vgw),&tmp); vgw.plusis(&tmp);
  vgw.normalize();
  // calculate angle in degrees:
  DL_Scalar inp=vgw.inprod(&vdw);
  if (inp>=1.0) ang=0; else
  if (inp<=-1.0) ang=180; else
                  ang=acos(inp)*(180/3.141592653);
  vgw.crossprod(&vdw,&tmp);
  if (tmp.inprod(&rdw)<0) ang=-ang;
  d_ang=ang-anglem; if (d_ang>=180) d_ang-=360; if (d_ang<-180) d_ang+=360;
  angle+=d_ang;
  anglem+=d_ang; if (anglem>=180) anglem-=360; if (anglem<-180) anglem+=360;
  return angle;
}
