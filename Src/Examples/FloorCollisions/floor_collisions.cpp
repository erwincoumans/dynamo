// A demo for the use of the Dynamo classes.
// It shows a chain of cubes that swings under influence of gravity.
// In this version the cubes are let go to collide with a floor. This
// shows the use and strengths and weaknesses of Dynamo's collision
// handling

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
		link->maxforce=25;
    }
    void Disconnect(){
      if (link) delete link;
      link=NULL;
    }

    // support for collision against a floor plane
    void DetectFloorCollision(double F/*floor height*/){
		if (companion->get_next_position()->y>F+sqrt(3)) return;
		// we're close to the floor; traverse all points to see which ones collide
		int i,j,k,nrcoll=0;
		for (i=-1;i<2;i+=2){
			for (j=-1;j<2;j+=2){
				for (k=-1;k<2;k+=2){
					DL_point pw,p(i,j,k);
					companion->new_toworld(&p,&pw);
					if (pw.y<F) { // p under the floor
						DL_vector v;
						companion->get_newvelocity(&p,&v);
						if (v.y<0) { // we have a collision
							pw.y=F;
							v.init(0,1,0);
							if (link && (i==-1) && (j==1) && (k==-1)) {
								// this point already collides through the point in
								// the cube above it
							}
							else {
								DL_collision *col=new DL_collision(companion,&p,NULL,&pw,&v,1);
								nrcoll++;
								if (nrcoll==3) i=j=k=2; // no need to look further
							}
						}
					}
				}
			}
		}
    }

    MyCube(DL_point& newpos){
    // constructor: its main task is creating and initialising the companion
      pos.assign(&newpos);
      orient.makeone();
      companion=new DL_dyna((void*)this);
      companion->set_mass(1);
      companion->set_inertiatensor(1,1,1);
      companion->set_velodamping(0.96); // introduce a bit of friction
	  companion->set_elasticity(0.8);
      link=NULL;
    }
    ~MyCube(){
      if (link) delete link;
      delete companion;
    }
};

MyCube* cube[NCUBES];

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
	virtual void do_collision_detection() {
		int i;
		for (i=0;i<NCUBES;i++) cube[i]->DetectFloorCollision((2.2-2*sqrt(3))*NCUBES-1);
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
  My_dyna_system_callbacks *dsc=new My_dyna_system_callbacks();
  DL_rungekutta2* my_int=new DL_rungekutta2();
  DL_dyna_system dsystem(dsc,my_int);
  DL_constraint_manager *constraints=new DL_constraint_manager();
  constraints->max_error=0.0001;

  DL_point pos(0.75*NCUBES,2.2*NCUBES,4*NCUBES);
  DL_vector vec(-2,-2,-2);

  for (i=0;i<NCUBES;i++){
    cube[i]=new MyCube(pos);
    pos.plusis(&vec);
  }
    
  cube[0]->ConnectTo(NULL);
  for (i=1;i<NCUBES;i++) cube[i]->ConnectTo(cube[i-1]);

  vec.init(0,-1,0);
  dsystem.set_gravity(&vec);
  my_int->set_stepsize(0.02);

  InitRender(dsystem,cube,NCUBES);

  // everything initialised. Now let the animation run:
#ifdef USINGDIRECTX
  while (RenderCubes(cube)) {
    dsystem.dynamics();
  }
#else
  RenderCubes();
#endif // USINGDIRECTX

  // all done: clean up:
  for (i=0;i<NCUBES;i++) delete cube[i];
  delete constraints; delete my_int; delete dsc;
}
