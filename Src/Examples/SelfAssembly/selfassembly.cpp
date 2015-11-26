// A demo for the use of the Dynamo classes.
// It shows a chain of cubes that swings under influence of gravity.
// Added to the basic example is the use of stiffness to allow slow
// self assembly of the chain.

// comment the following line to use OpenGL rendering
#define USINGDIRECTX

#include "rungekutta2.h"
#include "ptp.h"

// How many cubes make up the chain?:
#define NCUBES 15

// Here is the class for our cubes, these function mainly as an intermediate
// between the dyna companion, the control routine, and the rendering software.
class MyCube {
  public:
    // the attributes of the cube:
    DL_point pos;       // the position of the center of the cube
    DL_matrix orient;   // the orientation of the cube
    DL_dyna* companion; // the companion which does all the dynamics
                        // calculations for us
    DL_ptp* link;       // used to connect this cube to another

    // the functions needed for the DL_dyna_callback to copy position and
    // orientation to/from the companion:
    void get_new_geo_info(DL_geo *g){
      g->move(&pos,&orient);
    }
    void update_dyna_companion(DL_dyna *d){
      pos.assign(d->get_position());
      orient.assign(d->get_orientation());
    }

    // show an example of a constraint: this method uses a point-to-point
    // constraint to connect one of the cube's corners to a corner of the
    // other cube.
    void ConnectTo(MyCube *cube){
      DL_point posm(1,1,1),posp(-1,-1,-1);
      link=new DL_ptp();
      if (cube) link->init(companion,&posm,cube->companion,&posp);
      else { // connect to the current position of the corner in the world.
        DL_point posw;
        companion->to_world(&posm,&posw);
        link->init(companion,&posm,NULL,&posw);
      }
    }
    void Disconnect(){
      if (link) delete link;
      link=NULL;
    }

    MyCube(DL_point& newpos){
    // constructor: its main task is creating and initialising the companion
      pos.assign(&newpos);
      orient.makeone();
      companion=new DL_dyna((void*)this);
      companion->set_mass(1);
      companion->set_inertiatensor(1,1,1);
      companion->set_velodamping(0.995); // introduce a bit of friction
      link=NULL;
    }
    ~MyCube(){
      if (link) delete link;
      delete companion;
    }
};

// we need to implement the callbacks that allow the dyna system to commicate
// its calculated positions and orientations to us:
class My_dyna_system_callbacks : public DL_dyna_system_callbacks {
  public:
    virtual void get_new_geo_info(DL_geo *g){
      ((MyCube*)(g->get_companion()))->get_new_geo_info(g);
    }
    virtual void update_dyna_companion(DL_dyna *d){
      ((MyCube*)(d->get_companion()))->update_dyna_companion(d);
    }
};

// rendering is taken care of elsewhere:
#ifdef USINGDIRECTX
#include "render.cpp"					// use Microsoft Direct-X 
#else
#include "glrender.cpp"					// use OpenGL (Nagle)
#endif // USINGDIRECTX

// here comes the main control function: it initialises the dyna system
// and then creates a chain of cubes 
void main(){
  int i;
  MyCube* cube[NCUBES];
  My_dyna_system_callbacks *dsc=new My_dyna_system_callbacks();
  DL_m_integrator* my_int=new DL_rungekutta2();
  DL_dyna_system dsystem(dsc,my_int);
  DL_constraint_manager *constraints=new DL_constraint_manager();
  constraints->max_error=0.0001;

  DL_point pos(0.75*NCUBES,1.75*NCUBES,4*NCUBES);

  for (i=0;i<NCUBES;i++){
    cube[i]=new MyCube(pos);
  }
    
  cube[0]->ConnectTo(NULL);
  for (i=1;i<NCUBES;i++) cube[i]->ConnectTo(cube[i-1]);

  my_int->set_stepsize(0.02);

  InitRender(dsystem,cube,NCUBES);

  // everything initialised. Now start animating:

  // first set allow the chain to self assemble:
  for (i=0;i<NCUBES;i++) {
	  cube[i]->link->stiffness=0.005;
	  cube[i]->link->soft();
  }

#ifdef USINGDIRECTX
  while ((dsystem.time()<30) && RenderCubes(cube)) dsystem.dynamics();
#else
  // The OpenGL implementation will have to be changed here:
  RenderCubes();
#endif // USINGDIRECTX
  
  // now the chain has assembled: make the constraints hard again, and set
  // gravity to make the chain swing:

  for (i=0;i<NCUBES;i++) {
	  cube[i]->link->stiffness=1;
	  cube[i]->link->auto_softhard();
  }
  DL_vector grav(0,-1,0);
  dsystem.set_gravity(&grav);

  // swing along:

#ifdef USINGDIRECTX
  while (RenderCubes(cube)) {
    dsystem.dynamics();
    if (fabs(dsystem.time()-100)<0.01) cube[NCUBES/2]->Disconnect();
  }
#else
  RenderCubes();
#endif // USINGDIRECTX

  // all done: clean up:
  for (i=0;i<NCUBES;i++) delete cube[i];
  delete constraints; delete my_int; delete dsc;
}
