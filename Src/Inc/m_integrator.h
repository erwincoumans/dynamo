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
// filename: m_integrator.h
// description: class to derive ones own motion integrators from
// author: Bart Barenbrug   March '96
//

#ifndef DL_MINTEGRATORH
#define DL_MINTEGRATORH

#include "supvec.h"

class DL_dyna;

// ********************* //
// class DL_m_integrator //
// ********************* //

class DL_m_integrator {
  protected:
    DL_Scalar h;        // integration step
    DL_Scalar oldh;     //     ""  in the previous frame
    DL_Scalar halfh;    // h/2 
    DL_Scalar dfh;      // discretisation_factor*h (see dyna::ode() )
  public:
    /// for external (to DL) use:
    DL_Scalar stepsize();                     // get h;
    DL_Scalar halfstepsize();                 // get h/2;
    virtual void set_stepsize(DL_Scalar);     // set h;
    virtual DL_Scalar ast()=0; // get the avg time at which this integrator
                           // samples. Eg: since euler only samples at time t,
			   // it returns 0, but rungekutta4 samples at t,
			   // t+0.5*h and t+h, which averages to t+0.5*h.
			   // So rungekutta4 returns 0.5
    DL_m_integrator(){h=oldh=1.0; halfh=0.5; dfh=0.0;};   // constructor
    ~DL_m_integrator(){};                 // destructor

    /// for internal (DL) use only:
    virtual void integrate(DL_supvec*,DL_dyna*,DL_supvec*)=0;
                                          // do one integration step
    DL_Scalar old_stepsize(void){return oldh;}
    void shift_stepsize(void){oldh=h;}
};

#endif
