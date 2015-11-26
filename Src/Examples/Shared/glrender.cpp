//
//	glrender.cpp
//
//	OpenGL graphics for Dynamo dynamics system demo program
//
//	J. Nagle	May, 1999
//	nagle@animats.com
//	Released under the usual FSF General Public License.
//
//	This provides an alternative graphics library for
//	the DYNAMO dynamics system, for those who need to
//	use OpenGL rather than Microsoft's Direct-X.
//
//	Some minor changes have been made to dlexampl.cpp
//	to support this.
//
//	Related files: glrender.h, mycube.h
//
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>				
#include <GL/glaux.h>			
#include <process.h>
#include <vfw.h>
#include <new.h>
#include "ptp.h"
#include "matrix.h"
#include "algebra3.h"
//
//	Configuration
//
//
//	CLEAROVERLAYPLANES  -- clear overlay planes if enabled
//
//	This is a workaround for a bug in the Dynamic Pictures Oxygen driver for
//	NT 4.x systems, which doesn't clear the overlay planes when opening
//	a window on top of another window that uses overlay planes.
//
#define CLEAROVERLAYPLANES 1				// clear overlay planes if 1
//
class MyCube;
class DL_dyna_system;
//
//	CSimViewBackground  -- background info
//
class CSimViewBackground {
private:
	GLfloat fNear;									// near plane for fog/depth
	GLfloat fFar;									// far plane for fog/depth
	GLfloat fColor[4];
	bool	fHasFog;								// if fog or depth fading
	GLenum	fFogMode;								// GL fog mode
public:
	CSimViewBackground(double near, double far, const vec3& color, bool fog, bool depth);
	CSimViewBackground() { fHasFog = false; fColor[0]=fColor[1]=fColor[2]=0.0; fColor[3]=1.0;}
	void DrawBackground();							// build OpenGL background
};
//
//	CSimView  --  display class
//
//	This opens and updates an OpenGL window.
//
class CSimView
{
#if CLEAROVERLAYPLANES
	HGLRC	fOverlayContext;								// for clearing overlay plane
#endif // CLEAROVERLAYPLANES
	GLuint	fFontBase;										// base id for font
	CSimViewBackground	fBackground;						// background info
	vec3	fCameraInterest;								// camera params
	vec3	fCameraLocation;
	double	fCameraFieldOfView;
	//	DYNAMO objects
	int fCubeCount;											// number of cubes
	MyCube**	fCubes;										// the cubes
	DL_dyna_system* fSystem;								// system to simulate
private:
	void StartDraw();
	void ContentDraw();
	void FinishDraw();
	void DrawMsg(const char* s);							// output message
	void DrawLights();										// send lights to OpenGL
	void DrawLight(const int lightNum, const vec3& position, const vec3& interest, double coneAngle, vec3& color);
#if	CLEAROVERLAYPLANES
	void ClearOverlayPlanes();								// clear overlay planes
#endif // CLEAROVERLAYPLANES
public:
	static void DrawEvent();								// OpenGL Draw function
	static void CloseEvent();								// OpenGL has closed window
public:
	CSimView(DL_dyna_system& system, MyCube* cube[], int ncubes);												// use the factory
	virtual ~CSimView();										// destructor
	void GraphicsInit(const char window_title[]);			// initialization
	void Draw();											// OpenGL Draw function
	void SetCamera(const vec3& camera,  double zoom);		// set camera location
	void SetInterest(const vec3& camera);					// set camera interest
	void SetBackground(double near, double far, const vec3& color, bool fog, bool depth);
};

//
//	Configuration constants
//
//
//	Initial window size
const int kWindowWidth = 640;
const int kWindowHeight = 480;
//	Camera parameters
const GLdouble rollvector[] = {0.0,1.0,0.0 };			// "up" direction
const GLdouble kCameraNearPlane = 10.0;					// near plane
const GLdouble kCameraFarPlane = 65536*kCameraNearPlane;// far plane 
const vec3 kCamera(20,0,-20);							// camera location
const vec3 kInterest(10,0,60);							// camera looks here
const double kFieldOfView = 40.0;						// camera field of view
//	Lighting
const GLfloat kAmbientLight[4] = { 0.3f,0.3f,0.3f,1.0 };// ambient lighting
//	Depth cueing
const double kNearPlane = 10;
const double kFarPlane = 65536;
const vec3 kBackgroundColor(0,0,0.5f);					// dark blue
//	Font info (Win32 specific)
const char* kFontTypeface = "system";
const int kFontHeight = 20;
const DWORD kFontItalic = 0;
const int kFontWeight = 0;
//
void SetVertices(DL_point& pos, DL_matrix& orient);		// forward
//
//	The drawing object
//
CSimView* theDrawObject = NULL;	
//
//	
class MyVertex {
public:
	GLfloat x,y,z;
	GLint diffcolor;
};
MyVertex Vertex[24];

//
//	ToGLFloat3 - useful conversion
//
inline void ToGLfloat3(GLfloat out[3], const vec3& in)
{	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
}
//
//	ToGLFloat4 - useful conversion
//
inline void ToGLfloat4(GLfloat out[4], const vec3& in, double lastval)
{	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
	out[3] = lastval;			// the extra fourth parameter
}
//
//	GLCheckError -- check for openGL errors
//
static void GLCheckError ()
{
	GLenum err = glGetError ();		// Get error status
	if (err == GL_NO_ERROR) return;	// normal case
	fprintf (stderr, "OpenGL error: %s\n", gluErrorString (err));	// msg

	assert (err == GL_NO_ERROR);	// will trap 
	exit (1);						// fails 
}
//
//	Constructor
//
CSimView::CSimView(DL_dyna_system& system, MyCube* cube[], int ncubes)
{
	fSystem = &system;								// save link to dyna being simulated
	fCubes = cube;
	fCubeCount = ncubes;
#if CLEAROVERLAYPLANES
	fOverlayContext = NULL;							// no overlay context yet
#endif // CLEAROVERLAYPLANES
}
//
//	Destructor
//
CSimView::~CSimView()					
{

}

//
//	Draw  --  the usual Draw function
//
void CSimView::Draw()
{
	StartDraw();							// draw real info
	ContentDraw();
	FinishDraw();							// end drawing
	GLCheckError();
}

//
//	OpenGL support
//

//
//	glVertex  --  vertex for OpenGL using vec3.
//
static void glVertex(const vec3& pt)
{	glVertex3f(pt[0],pt[1],pt[2]);				// do it
}
//
//	DrawLine  --  draw a line between two 3D points.
//
//	Immediate.
//
void DrawLine(const vec3& pt1, const vec3& pt2)
{
	glPushMatrix();
	glLoadIdentity();
	glDisable (GL_LIGHTING);
	glLineWidth(2.0);
	glColor3ub(200, 20, 20);
	glBegin(GL_LINES);
	glVertex(pt1);
	glVertex(pt2);
	glEnd();
	glEnable (GL_LIGHTING);
	glPopMatrix();
}
//
//	ReshapeViewport  --  resizes the window
//
//	All we do here is resize the viewport.  The matrices are computed in
//	StartDraw.
//
static void CALLBACK 
ReshapeViewport (int width, int height)
{
	if (height < 1) height = 1;			// OpenGL doesn't like zero-size viewports
	if (width < 1) width = 1;			// although NT lets you resize to zero.
	glViewport (0, 0, width, height);
}
//
const WORD Idx[36]=
{ 0,  1,  3,
  1,  2,  3,
  4,  7,  5,
  5,  7,  6,
  9,  8, 10,
  8, 11, 10,
 12, 13, 15,
 14, 15, 13,
 18, 17, 16,
 16, 19, 18,
 23, 20, 21,
 21, 22, 23
};
//
//	DrawTriangle  --  draw one triangle with OpenGL.  
//
static void DrawTriangle(const MyVertex& v0, const MyVertex& v1, const MyVertex& v2)
{
	//	Calculate edge vectors
	GLfloat edge0[3], edge1[3];
	edge0[0] = v0.x - v1.x;
	edge0[1] = v0.y - v1.y;
	edge0[2] = v0.z - v1.z;
	edge1[0] = v1.x - v2.x;
	edge1[1] = v1.y - v2.y;
	edge1[2] = v1.z - v2.z;
	//	Compute the normal, which is the cross product of the edge vectors
	GLfloat norm[3];
	norm[0] = edge0[1]*edge1[2] - edge0[2]*edge1[1];
	norm[1] = edge0[2]*edge1[0] - edge0[0]*edge1[2];
	norm[2] = edge0[0]*edge1[1] - edge0[1]*edge1[0];
	//	Normalize the normal
	const double length = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
	norm[0] /= length;
	norm[1] /= length;
	norm[2] /= length;
	glNormal3f(norm[0],norm[1],norm[2]);		// specify the normal

	//	Draw the vertices
	glVertex3f(v0.x,v0.y,v0.z);		// draw it
	glVertex3f(v1.x,v1.y,v1.z);		// draw it
	glVertex3f(v2.x,v2.y,v2.z);		// draw it
}
//
//	ContentDraw
//
//	Enter with ModelView and Projection matrices set.
//
void CSimView::ContentDraw()
{
	//	Run the simulation
	static int framenumber;
	const int kShowInterval = 10;			// show every 10th step
	for (int step=0; step < kShowInterval; step++)
	{	fSystem->dynamics();
		framenumber++;
	}
	//	Draw all active polytopes.
	vec4 color(.5,.5,.5,1);						// white
	GLfloat params[4];
	params[0] = color[0];
	params[1] = color[1];
	params[2] = color[2];
	params[3] = color[3];
	glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,params);	// set material property
	////glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,params);	// set material property
	glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,params);	// set material property
	for (int i=0;i<fCubeCount;i++)
	{	SetVertices(fCubes[i]->pos, fCubes[i]->orient);
		glBegin(GL_TRIANGLES);
		for (int j=0; j<12; j++)				// draw the triangles of the cube
		{	const int ix = j*3;					// do a face
			////DrawTriangle(Vertex[Idx[ix]],Vertex[Idx[ix+1]],Vertex[Idx[ix+2]]);
			DrawTriangle(Vertex[Idx[ix+2]],Vertex[Idx[ix+1]],Vertex[Idx[ix]]);

		}
		glEnd();
	}
	char ss[50];									// progress indicator
	_snprintf(ss,sizeof(ss),"Frame %d",framenumber);		// edit
	DrawMsg(ss);									// display

}
//
//	StartDraw()  -- start drawing
//
void CSimView::StartDraw()
{
	const float color = float(0.3);
	//	Compute projection matrix
	glMatrixMode(GL_PROJECTION);					// work on model/view xform
	GLint viewportParams[4];						// x,y,w,h
	glGetIntegerv(GL_VIEWPORT,viewportParams);		// get viewpoint info
	GLfloat width = viewportParams[2];				// window width
	GLfloat height = viewportParams[3];				// window height
	if (height < 1.0) height = 1.0;					// avoid divide by zero
	GLfloat aspectRatio = width / height;			// compute aspect
	glLoadIdentity ();
	{	
		gluPerspective(fCameraFieldOfView,aspectRatio,kCameraNearPlane,kCameraFarPlane);	// set perspective
		//	Computer modelview matrix
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(fCameraLocation[0],fCameraLocation[1],fCameraLocation[2],
			fCameraInterest[0],fCameraInterest[1],fCameraInterest[2],
			rollvector[0],rollvector[1],rollvector[2]);
	}
	fBackground.DrawBackground();					// clear to background color, set fog
#if CLEAROVERLAYPLANES
	ClearOverlayPlanes();							// clear overlay planes
#endif // CLEAROVERLAYPLANES
	DrawLights();									// put in all the lights
	//	Transparency
	////glEnable(GL_BLEND);
	////glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
}
//
//	FinishDraw()  -- done with drawing
//
void CSimView::FinishDraw()
{	glFlush();
	auxSwapBuffers ();
	GLCheckError ();				// OpenGL happy?
}
//
//	DrawMsg  -- draw a text message
//
//	Chars must be in the printable range 32-127
//
void CSimView::DrawMsg(const char* msg)
{	assert(msg);
	if (!fFontBase) return;
	glPushAttrib(GL_LIST_BIT);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);					// work on projection matrix
	glPushMatrix();									// save view transform
	glDisable (GL_LIGHTING);						// no lighting
		glLoadIdentity();							// directly to screen
		gluOrtho2D(0,639,0,480);					// ***TEMP TEST***
 		glColor3f(1, 0, 0);							// write in red
		glListBase(fFontBase-32);					// ASCII offset
		glRasterPos2i(kFontHeight,kFontHeight);		// One char height from bottom, equallly from left
		glCallLists(strlen(msg),GL_UNSIGNED_BYTE,msg); // draw text
		GLCheckError();
	glEnable (GL_LIGHTING);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glPopAttrib();
	GLCheckError();
}
//
//	GLDrawEvent  --  OpenGL calls this to force a draw
//
static void CALLBACK GLDrawEvent(void)
{
	if (theDrawObject) theDrawObject->Draw();		// do it
}

//
//	GraphicsEventLoop
//
void 
GraphicsEventLoop (void (CALLBACK* eventFn) (void))	
{	GLCheckError();
	auxReshapeFunc(ReshapeViewport);// identify to aux library
	GLCheckError();
	auxIdleFunc (eventFn);			// no idle processing
	GLCheckError();	
	auxMainLoop (eventFn);			// all the work is done here
	GLenum err = glGetError ();		// Get error status 
	if (err == GL_NO_ERROR) return;	// normal
	if (err == GL_INVALID_OPERATION) 
	{	fprintf(stderr,"Stopping because the window was closed.\n");
		return;
	}
	GLCheckError();
}
//
//	GetScreenSize  --  get screen size from system
//
//	System-dependent
//
static void GetScreenSize(GLint& x, GLint& y)
{
	x = GetSystemMetrics(SM_CXSCREEN);					// screen X size, pixels	
	y = GetSystemMetrics(SM_CYSCREEN);					// screen Y size, pixels
}
#if CLEAROVERLAYPLANES
//
//	ClearOverlayPlanes  --  clears overlay planes, if any.
//
//	A workaround.  Some OpenGL drivers are too dumb to do this automatically.
//
void CSimView::ClearOverlayPlanes()
{
	HDC hdc = wglGetCurrentDC();					// get current device context
	if (!hdc) return;								// no context, fails
	HGLRC oldContext = wglGetCurrentContext();		// get current OpenGL context, if any
	if (!oldContext) return;						// no current OpenGL context, fails
	if (!fOverlayContext) return;					// no overlay plane available
	 
	// make it the calling thread's current rendering context 
	BOOL stat = wglMakeCurrent(hdc, fOverlayContext); 
	if (!stat) return;
 
	//	Clear the overlay buffer.
	glClearColor(1,1,1,0.0);						// clear to white
 	glClear(GL_COLOR_BUFFER_BIT);					// clear buffer
	// Done with the overlay plane rendering context. 
	wglMakeCurrent(hdc, oldContext);				// back to old OpenGL context
}
#endif // CLEAROVERLAYPLANES
//
//	GraphicsInit  --  open an OpenGL window
//
void CSimView::GraphicsInit(const char window_title[])
{	
	auxInitDisplayMode (AUX_DOUBLE | AUX_RGBA | AUX_DEPTH);
	GLint maxwidth,maxheight;
	GetScreenSize(maxwidth,maxheight);					// get screen size info
	GLint left = (maxwidth - kWindowWidth)/2;			// center on screen
	GLint top = (maxheight - kWindowHeight)/2;
	assert(left >= 0);									// must fit on screen
	assert(top >= 0);
	auxInitPosition (left, top, kWindowWidth, kWindowHeight);
	glViewport(0, 0, kWindowWidth, kWindowHeight);
	auxInitWindow (window_title);
	glLineWidth (2.0);
	glEnable (GL_LIGHTING);
	////glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,0.0);				// local viewer lighting mode
	glDepthFunc(GL_LEQUAL);
	glEnable (GL_DEPTH_TEST);
	glShadeModel (GL_SMOOTH);
	//	Font initialization  -- WIN32 specific
	fFontBase = glGenLists(96);						// 96 printable ASCII chars
	if (fFontBase)									// per OpenGL SuperBible
	{	HFONT font = CreateFont(kFontHeight, 0, 0, 0, kFontWeight, kFontItalic, FALSE, FALSE,
		ANSI_CHARSET,OUT_TT_PRECIS,
		CLIP_DEFAULT_PRECIS,DRAFT_QUALITY,
		FIXED_PITCH,kFontTypeface);
		if (font)
		{	HDC hdc = wglGetCurrentDC();			// get current device context
			SelectObject(hdc,font);
			wglUseFontBitmaps(hdc,32,96,fFontBase);
		}
	}
	//	Clearing
	glClearColor(0.0,0.0,0.0,0);					// clear to black
	glClear(GL_COLOR_BUFFER_BIT);					// clear
#if CLEAROVERLAYPLANES								// if clearing overlay planes
	HDC hdc = wglGetCurrentDC();					// get current device context
	if (hdc)
	{	fOverlayContext = wglCreateLayerContext(hdc,PFD_OVERLAY_PLANE);	 // save plane handle
		ClearOverlayPlanes();						// clear the overlay planes, if any
	}
#endif // CLEAROVERLAYPLANES
}

//
//	SetCamera -- set camera location and field of view
//
void CSimView::SetCamera(const vec3& camera,  double fov)
{	
	fCameraFieldOfView = fov;						// field of view (degrees)
	fCameraLocation = camera;
}
//
//	SetInterest --  look at this point
//
void CSimView::SetInterest(const vec3& interest)
{
	fCameraInterest = interest;
}
//
//	SetBackground  --  set fog and depth-fading info
//
void CSimView::SetBackground(double nearplane, double farplane, const vec3& color, bool fog, bool depth)
{
	fBackground = CSimViewBackground(nearplane,farplane,color,fog,depth);	// save
}
//
//	DrawLight  --  draw a single light
//
//	Point light:			cone angle = 360
//	Spotlight:				cone angle = cone angle
//	Infinite light or sun:	cone angle = 0
//
void CSimView::DrawLight(const int lightNum, const vec3& positionin, const vec3& interest, double coneAngle, vec3& color)
{
	const double kBrightness = 1.0;								// overall brightness ***TEMP***
	const double kPointLightBrightness = 0.1;					// point light brightness ***TEMP***
	bool pointsource = coneAngle > 0.0001;						// true if directional light
	double spotCutoff = coneAngle;								// open GL spot cutoff
	if (spotCutoff < 0.001 || spotCutoff >= 90.0)				// if no cone angle
	{	spotCutoff = 180;										// turn off cone
	}
	GLfloat position[4];
	ToGLfloat4(position,positionin,pointsource ? 1.0 : 0.0);
	vec3 dir = interest-positionin;	// vector for light direction
	double len = dir.length();		// prepare to normalize
	if (len <= 0.000001)
	{	dir = vec3(0,1,0);			// light from above
	} else {
		dir /= len;					// normalize dir
	}
	vec3 dimcolor(color*kBrightness);	// adjust brightness
	GLfloat spotDirection[3];		// direction for spotlight
	ToGLfloat3(spotDirection,dir);
	GLfloat wcolor[4];
	ToGLfloat4(wcolor,color,1.0);	// convert color
	glLightfv(GL_LIGHT0+lightNum,GL_POSITION,position);		// position
	glLightfv(GL_LIGHT0+lightNum,GL_SPOT_DIRECTION,spotDirection);	// direction
	glLightf(GL_LIGHT0+lightNum,GL_SPOT_CUTOFF,spotCutoff);	// spot coverage
															// NO AMBIENT PART
	glLightfv(GL_LIGHT0+lightNum,GL_DIFFUSE,wcolor);		// diffuse light part
	glLightfv(GL_LIGHT0+lightNum,GL_SPECULAR,wcolor);		// specular light part
	glEnable(GL_LIGHT0+lightNum);							// enable this light
	GLCheckError();											// detect trouble
}
//
//	DrawLights  --  generate the lights in OpenGL
//
void CSimView::DrawLights()
{	//	Generate all the lights, up to the OpenGL limit
	GLCheckError();
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,kAmbientLight);	// ambient light
	DrawLight(0,vec3(500,500,-500),vec3(0,0,0),0,vec3(1,1,1));
	DrawLight(1,vec3(-500,500,-500),vec3(0,0,0),0,vec3(1,1,1));
	GLCheckError();
}
//
//	Class CSimViewBackground
//
CSimViewBackground::CSimViewBackground(double nearlim, double farlim, const vec3& color, bool fog, bool depth)
{	ToGLfloat4(fColor,color,1.0);			// convert color
	fNear = nearlim;						// near and far planes
	fFar = farlim;
	fHasFog = false;
	if (depth)								// depth fading is linear fog
	{	fFogMode = GL_LINEAR;
		fHasFog = true;
	}
	if (fog)								// fog is exponential fog
	{	fFogMode = GL_EXP;
		fHasFog = true;
	}
}
//
//	DrawBackground -- set all background light params
//
void CSimViewBackground::DrawBackground()
{
	glClearColor(fColor[0],fColor[1],fColor[2],1.0);	// background color
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);// clear all to it
	if (fHasFog)										// if fog or depth cueing
	{	glEnable(GL_FOG);								// set fog params
		glFogf(GL_FOG_MODE,fFogMode);
		glFogfv(GL_FOG_COLOR,fColor);
		glFogf(GL_FOG_START,fNear);
		glFogf(GL_FOG_END,fFar);
	}
}
//
//	Interfaces with the DYNAMO example program
//
//
//	SetVertices  --  position and orient cube at indicated location.
//
static void SetVertices(DL_point& pos, DL_matrix& orient){
	////DL_matrix orient;
	////orientin.transpose(&orient);
	Vertex[16].x=Vertex[12].x=Vertex[0].x=(float)(pos.x-orient.c0.x-orient.c1.x+orient.c2.x);
	Vertex[16].y=Vertex[12].y=Vertex[0].y=(float)(pos.y-orient.c0.y-orient.c1.y+orient.c2.y);
	Vertex[16].z=Vertex[12].z=Vertex[0].z=(float)(pos.z-orient.c0.z-orient.c1.z+orient.c2.z);
	
	Vertex[17].x=Vertex[8].x=Vertex[1].x=(float)(pos.x+orient.c0.x-orient.c1.x+orient.c2.x);
	Vertex[17].y=Vertex[8].y=Vertex[1].y=(float)(pos.y+orient.c0.y-orient.c1.y+orient.c2.y);
	Vertex[17].z=Vertex[8].z=Vertex[1].z=(float)(pos.z+orient.c0.z-orient.c1.z+orient.c2.z);
	
	Vertex[21].x=Vertex[9].x=Vertex[2].x=(float)(pos.x+orient.c0.x-orient.c1.x-orient.c2.x);
	Vertex[21].y=Vertex[9].y=Vertex[2].y=(float)(pos.y+orient.c0.y-orient.c1.y-orient.c2.y);
	Vertex[21].z=Vertex[9].z=Vertex[2].z=(float)(pos.z+orient.c0.z-orient.c1.z-orient.c2.z);
	
	Vertex[20].x=Vertex[13].x=Vertex[3].x=(float)(pos.x-orient.c0.x-orient.c1.x-orient.c2.x);
	Vertex[20].y=Vertex[13].y=Vertex[3].y=(float)(pos.y-orient.c0.y-orient.c1.y-orient.c2.y);
	Vertex[20].z=Vertex[13].z=Vertex[3].z=(float)(pos.z-orient.c0.z-orient.c1.z-orient.c2.z);
	
	Vertex[19].x=Vertex[15].x=Vertex[4].x=(float)(pos.x-orient.c0.x+orient.c1.x+orient.c2.x);
	Vertex[19].y=Vertex[15].y=Vertex[4].y=(float)(pos.y-orient.c0.y+orient.c1.y+orient.c2.y);
	Vertex[19].z=Vertex[15].z=Vertex[4].z=(float)(pos.z-orient.c0.z+orient.c1.z+orient.c2.z);
	
	Vertex[18].x=Vertex[11].x=Vertex[5].x=(float)(pos.x+orient.c0.x+orient.c1.x+orient.c2.x);
	Vertex[18].y=Vertex[11].y=Vertex[5].y=(float)(pos.y+orient.c0.y+orient.c1.y+orient.c2.y);
	Vertex[18].z=Vertex[11].z=Vertex[5].z=(float)(pos.z+orient.c0.z+orient.c1.z+orient.c2.z);
	
	Vertex[22].x=Vertex[10].x=Vertex[6].x=(float)(pos.x+orient.c0.x+orient.c1.x-orient.c2.x);
	Vertex[22].y=Vertex[10].y=Vertex[6].y=(float)(pos.y+orient.c0.y+orient.c1.y-orient.c2.y);
	Vertex[22].z=Vertex[10].z=Vertex[6].z=(float)(pos.z+orient.c0.z+orient.c1.z-orient.c2.z);
	
	Vertex[23].x=Vertex[14].x=Vertex[7].x=(float)(pos.x-orient.c0.x+orient.c1.x-orient.c2.x);
	Vertex[23].y=Vertex[14].y=Vertex[7].y=(float)(pos.y-orient.c0.y+orient.c1.y-orient.c2.y);
	Vertex[23].z=Vertex[14].z=Vertex[7].z=(float)(pos.z-orient.c0.z+orient.c1.z-orient.c2.z);
}


//
//	RenderCubes -- Called from the example
//
//	Called once, and does all the work.
//
void RenderCubes()
{
	GraphicsEventLoop(GLDrawEvent);
}
//
//	InitRender --  initialize OpenGL rendering
//
int InitRender(DL_dyna_system& system, MyCube* cube[], int ncubes)
{
	assert(theDrawObject == NULL);						// no draw object yet
	theDrawObject = new CSimView(system,cube,ncubes);	// create a draw object
	theDrawObject->GraphicsInit("Dynamo example");		// open the window
	theDrawObject->SetInterest(kInterest);				// look here
	theDrawObject->SetCamera(kCamera,kFieldOfView);		// camera here
	theDrawObject->SetBackground(kNearPlane,kFarPlane,kBackgroundColor,false,true);
	return(0);
}
//
//	Dummy WinMain, for compatibility with dlexampl.cpp
//
void main();											// external
int
WINAPI
WinMain(
    HINSTANCE hInstance,
    HINSTANCE hPrevInstance,
    LPSTR lpCmdLine,
    int nShowCmd
    )
{	main();
	return(0);
}



