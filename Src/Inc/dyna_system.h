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

// filename: dyna_system.h

// description: system class controlling dynamic motion calculations

#ifndef DL_DYNASYSTEMH
#define DL_DYNASYSTEMH

#include <stdarg.h>
#include <stdio.h>
#include "m_integrator.h"
#include "geo.h"
#include "controller.h"
#include "list.h"
#include "force_drawer.h"

class DL_dyna;

// ****************************** //
// class DL_dyna_system_callbacks //
// ****************************** //

class DL_dyna_system_callbacks : public DL_force_drawer {
  public:
    virtual void get_first_geo_info(DL_geo *g){
			// to be (optionally) overriden by a descendant:
			// implements the initial copying of the
			// position and velocity of the geo's companion
			// to the DL_geo
	// default implementation:
	// get the position and orientation from g:
	get_new_geo_info(g);
	// copy this new location to the old location:
	g->set_position(g->DL_geo::get_next_position());
	g->set_orientation(g->DL_geo::get_next_orientation());
	// and set the velocities to zero
	DL_vector v(0,0,0);
	g->set_velocity(&v);	g->set_next_velocity(&v);
	g->set_angvelocity(&v); g->set_next_angvelocity(&v);
    };
    virtual void get_new_geo_info(DL_geo*)=0;
                        // to be overriden by a descendant:
		        // implements copying the new
			// position of the geo's companion
			// to the DL_geo (using DL_geo::move(..)
			// for example)
    virtual void update_dyna_companion(DL_dyna*)=0;
                        // to be overridden by a descendent:
			// that implements copying the dyna's
			// new position and orientation to
			// the dyna's companion
    virtual void check_inertiatensor(DL_dyna*){};
                        // to be overridden by a descendent:
			// implements making sure the inertia
			// tensor of the DL_dyna is ok
    virtual void do_collision_detection(void){};
                        // do the collsion detection, which -upon detection
			// of a collision- creates and initialises a collision
			// constraint which will take care of the collision
			// handling
			// by default no collision detection is called
    virtual void Msg(char *fmt, ...){
      va_list args;
  
      va_start(args,fmt);
      vfprintf(stderr,fmt,args);
      va_end(args);
    }

};

// ******************** //
// class DL_dyna_system //
// ******************** //

class DL_dyna_system {
  protected:
    int frame_nr;                // the frame number
    DL_Scalar curtime;           // the current time
    DL_vector gravity;           // global gravity
    DL_m_integrator *integrator; // the motion integrator used by all dyna's
    DL_dyna_system_callbacks *companion;
                                 // the companion object for interfacing with
				 // the outside user
    DL_List geos;                // these are the geo's that are
		 	         // of interest to the dyna_system.
    DL_List controllers;         // these are the controllers that are
				 // managed by the dyna_system.
    boolean show_con_forces;     // show the controller forces or not
  public:
    /// for external (to DL) use:
    DL_dyna_system_callbacks* get_companion(void){ return companion; };

    void dynamics(void);                // do the dynamics (entry point)
    
    DL_geo* register_geo(void*);        // make sure there is a geo with the
                                        // supplied companion in the geos-list
    void remove_geo(DL_geo*);           // remove the geo from the geos.
    
    void set_gravity(DL_vector*);       // set gravity
    void get_gravity(DL_vector*);       // get gravity;
    void set_integrator(DL_m_integrator*); // set int
    DL_m_integrator* get_integrator(void); // get int
    DL_Scalar kinenergy();              // return sum of kinetic
                                        // energy of all dyna's
    DL_Scalar potenergy();              // return sum of potential
                                        // energy of all dyna's
    DL_Scalar totenergy();              // return sum of total
                                        // energy of all dyna's

    int frame_number() { return frame_nr; };
    DL_Scalar time() { return curtime; };

    void show_controller_forces();
    void hide_controller_forces();
    boolean showing_controller_forces(){return show_con_forces;};

    DL_dyna_system(DL_dyna_system_callbacks*,DL_m_integrator*);
                       // constructor
    //~DL_dyna_system(){}; // destructor

    /// for internal (DL) use only:
    void register_dyna(DL_dyna*);       // add the DL_dyna to the dynas-list
    void remove_dyna(DL_dyna*);         // remove the dyna from the dynas.
    void add_controller(DL_controller*); // add this controller to the list
    void rem_controller(DL_controller*); // remove this controller from the list
    void update_dyna_companions();       // update positions/orientations etc. for all dynas
    
    DL_Scalar newkinenergy();           // return sum of new kinetic
                                        // energy of all dyna's
    DL_Scalar newpotenergy();           // return sum of new potential
                                        // energy of all dyna's
    DL_Scalar newtotenergy();           // return sum of new total
                                        // energy of all dyna's


    // the next attribute is *temporarily* public for the makeshift collision
    // detector in the GDP...
    DL_List dynas;                      // these are the dyna's that are
				        // managed by this dyna_system.
};
				       
extern DL_dyna_system* DL_dsystem;

#endif
