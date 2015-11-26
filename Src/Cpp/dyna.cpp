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
// filename     : dyna.cpp
// description	: non-inline methods of class DL_dyna
// author	: Bart Barenbrug   March 1996
//

#include "dyna.h"
#include "dyna_system.h"

//#define DEBUG
//#define DEBUG2

// ************************** //
// non-inline member fuctions //
// ************************** //

void DL_dyna::assign(DL_dyna* d, void *comp) {
  if (d) {
    DL_geo::assign(d,comp);
    oneD=d->oneD;
    J.assign(&(d->J)); 
    Jinv.assign(&(d->Jinv)); 
    F.assign(&(d->F));
    M.assign(&(d->M));
    // assign the forces list:
    {
      forces.delete_all();
      DL_Mpair *tmp=(DL_Mpair*)d->forces.getfirst();
      while (tmp) {
        DL_Mpair *tmp2=new DL_Mpair;
	tmp2->force.assign(&(tmp->force));
	tmp2->rho.assign(&(tmp->rho));
	forces.addelem(tmp2);
	tmp=(DL_Mpair*)(d->forces.getnext(tmp));
      }
    }
    nextmstate.assign(&(d->nextmstate));
    mstateimp.assign(&(d->mstateimp));
    Fuptodate=d->Fuptodate;
    Muptodate=d->Muptodate;
    velodamping=d->velodamping;
    totalmass=d->totalmass;
    totalmass_inv=d->totalmass_inv;
  }
}

void DL_dyna::prepare_for_next_frame(void) {
  // integrate the motion state to obtain nextmstate
  integrate();  

  // shift statevectors and reset total applied force/torque
  nextmstate.q.normalize();
  nextmstate.q2A(this);
  mstate.assign(&nextmstate);
  Fuptodate=Muptodate=FALSE;
  F.init(0,0,0); M.init(0,0,0);
  forces.delete_all();
  matrixcache1empty=matrixcache2empty=TRUE;

  // apply velocity damping.
  DL_Scalar vdh=pow(velodamping,DL_dsystem->get_integrator()->stepsize());
		    // (account for stepsize);
  mstate.v.timesis(vdh);
  mstate.w.timesis(vdh);
  mstateimp.assign(&mstate);

#ifdef DEBUG
 extern int lookstime;
  dystem->get_companion()->Msg("[%d] angular velocity of %s: %f %f %f, length: %f\n",
       lookstime, getname(), mstate.w.x,  mstate.w.y, mstate.w.z, mstate.w.norm() );
  dystem->get_companion()->Msg("center of mass of %s: %f %f %f\n", getname(), mstate.z.x, mstate.z.y, mstate.z.z );
  dystem->get_companion()->Msg("orientation of %s:\n %f %f %f\n %f %f %f\n %f %f %f\n",
      getname(),
      mstate.A.c0->x, mstate.A.c1->x, mstate.A.c2->x,
      mstate.A.c0->y, mstate.A.c1->y, mstate.A.c2->y,
      mstate.A.c0->z, mstate.A.c1->z, mstate.A.c2->z
     );
#endif
}

void DL_dyna::begintest(void) {
  if (testing) return; // already copied everything
  testing=TRUE;
  FSave.assign(&F);
  MSave.assign(&M);
  FuptodateSave=Fuptodate;
  MuptodateSave=Muptodate;
  mstateimpSave.assign(&mstateimp);
  if (Fuptodate || Muptodate) {
    mstateSave.assign(&nextmstate);
    // save the forces list:
    {
      forcesSave.delete_all();
      DL_Mpair *tmp=(DL_Mpair*)forces.getfirst();
      while (tmp) {
        DL_Mpair *tmp2=new DL_Mpair;
	tmp2->force.assign(&(tmp->force));
	tmp2->rho.assign(&(tmp->rho));
	forcesSave.addelem(tmp2);
	tmp=(DL_Mpair*)forces.getnext(tmp);
      }
    }    
  }
}

void DL_dyna::endtest(void) {
  if (!testing) return; // already copied everything back
  testing=FALSE;
  F.assign(&FSave);
  M.assign(&MSave);
  Fuptodate=FuptodateSave;
  Muptodate=MuptodateSave;
  mstateimp.assign(&mstateimpSave);
  if (Fuptodate || Muptodate) {
    nextmstate.assign(&mstateSave);
    // restore the forces list:
    {
      forces.delete_all();
      DL_Mpair *tmp=(DL_Mpair*)forcesSave.getfirst();
      while (tmp) {
        DL_Mpair *tmp2=new DL_Mpair;
	tmp2->force.assign(&(tmp->force));
	tmp2->rho.assign(&(tmp->rho));
	forces.addelem(tmp2);
	tmp=(DL_Mpair*)forcesSave.getnext(tmp);
      }
    }    
  }
}

DL_Scalar DL_dyna::kinenergy(void) {
  DL_vector vtmp;
  vtmp.x=mstate.A.c0.inprod(&(mstate.w));
  vtmp.y=mstate.A.c1.inprod(&(mstate.w));
  vtmp.z=mstate.A.c2.inprod(&(mstate.w));
  // vtmp=A^T*omega
  return 0.5*(totalmass*mstate.v.inprod(&(mstate.v))+
              vtmp.x*vtmp.x*J.x+
	      vtmp.y*vtmp.y*J.y+
	      vtmp.z*vtmp.z*J.z);
}

DL_Scalar DL_dyna::potenergy(void) {
  DL_vector zv;
  mstate.z.tovector(&zv);
  return -Fexternal.inprod(&zv);
}

DL_Scalar DL_dyna::totenergy(void) {
  return kinenergy()+potenergy();
}

DL_Scalar DL_dyna::impkinenergy(void) {
  DL_vector vtmp;
  vtmp.x=mstateimp.A.c0.inprod(&(mstateimp.w));
  vtmp.y=mstateimp.A.c1.inprod(&(mstateimp.w));
  vtmp.z=mstateimp.A.c2.inprod(&(mstateimp.w));
  // vtmp=A^T*omega
  return 0.5*(mstateimp.v.inprod(&(mstateimp.v))*totalmass+
              vtmp.x*vtmp.x*J.x+
	      vtmp.y*vtmp.y*J.y+
	      vtmp.z*vtmp.z*J.z);
}

DL_Scalar DL_dyna::imppotenergy(void) {
  DL_vector zv;
  mstateimp.z.tovector(&zv);
  return -Fexternal.inprod(&zv);
}

DL_Scalar DL_dyna::imptotenergy(void) {
  return impkinenergy()+imppotenergy();
}

DL_Scalar DL_dyna::newkinenergy(void) {
  DL_vector vtmp;
  integrate();
  vtmp.x=nextmstate.A.c0.inprod(&(nextmstate.w));
  vtmp.y=nextmstate.A.c1.inprod(&(nextmstate.w));
  vtmp.z=nextmstate.A.c2.inprod(&(nextmstate.w));
  // vtmp=A^T*omega
  return 0.5*(totalmass*nextmstate.v.inprod(&(nextmstate.v))+
              vtmp.x*vtmp.x*J.x+
	      vtmp.y*vtmp.y*J.y+
	      vtmp.z*vtmp.z*J.z);
}

DL_Scalar DL_dyna::newpotenergy(void) {
  DL_vector zv;
  integrate();
  nextmstate.z.tovector(&zv);
  return -Fexternal.inprod(&zv);
}

DL_Scalar DL_dyna::newtotenergy(void) {
  return newkinenergy()+newpotenergy();
}

#undef DEBUG
