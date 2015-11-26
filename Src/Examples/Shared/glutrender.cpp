// render functionality using Glut, donated by Vincent Harvey

#include <GL/glut.h>
#include "geo.h"

static DL_dyna_system* Rsystem;
static MyCube** Rcube;
static int Rncubes;

void getGLMatrix (MyCube* c, float matrix[16])
{
    matrix[0] = c->orient.c0.x;
    matrix[1] = c->orient.c0.y;
    matrix[2] = c->orient.c0.z;
    matrix[3] = 0.0;

    matrix[4] = c->orient.c1.x;
    matrix[5] = c->orient.c1.y;
    matrix[6] = c->orient.c1.z;
    matrix[7] = 0.0;

    matrix[8] = c->orient.c2.x;
    matrix[9] = c->orient.c2.y;
    matrix[10] = c->orient.c2.z;
    matrix[11] = 0.0;

    matrix[12] = c->pos.x;
    matrix[13] = c->pos.y;
    matrix[14] = c->pos.z;
    matrix[15] = 1.0;
}

void draw ()
{
    float matrix[16];
    
    Rsystem->dynamics();
    
    glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glLoadIdentity ();

    GLfloat light_pos[] = { 0, 50, -50, 1 };
    GLfloat light_ambient[] = { 0.5, 0.5, 1.0, 1.0 };

    glLightfv ( GL_LIGHT0, GL_AMBIENT, light_ambient );
    glLightfv ( GL_LIGHT0, GL_POSITION, light_pos );
    glColor3f ( 0.5, 0.8, 0.9 );
    
    glColor3f ( 1.0, 1.0, 1.0 );
    
    // This should get the view right in the center of the action,
    // not to close, to to far, depending how many cubes there are.
    glTranslatef (-1.0 * Rncubes, 0.0, -10 * Rncubes);

    for (int i = 0; i < Rncubes; i++)
    {
        glPushMatrix ();
        getGLMatrix (Rcube[i],matrix);
        glMultMatrixf (matrix);
        glutSolidCube (2);
        glPopMatrix ();
    }

    glutSwapBuffers ();
}

void reshapeWindow ( int width, int height )
{
    glViewport ( 0, 0, width, height );
    glMatrixMode ( GL_PROJECTION );
    glLoadIdentity ();
    glFrustum ( -1.0, 1.0, -1.0, 1.0, 2.0, 1000.0 );
    glMatrixMode ( GL_MODELVIEW );
    glLoadIdentity ();
}

void idle ()
{
    glutPostRedisplay ();
}

int InitRender (DL_dyna_system& Isystem, MyCube* Icube[], int Incubes)
{
    Rsystem = &Isystem;
    Rcube = Icube;
    Rncubes = Incubes;
    
    glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
    glutInitWindowSize ( 450, 450 );
    glutInitWindowPosition ( 0, 0 );
    glutCreateWindow ( "Dynamo Animation Viewer" );

    glClearColor ( 0.0, 0.0, 0.0, 0.0 );
    glEnable ( GL_LIGHTING );
    glEnable ( GL_LIGHT0 );
    glEnable ( GL_DEPTH_TEST );
    glDepthFunc ( GL_LEQUAL );
    glShadeModel ( GL_SMOOTH );

    glutReshapeFunc ( reshapeWindow );
    glutIdleFunc ( idle );
    glutDisplayFunc ( draw );

    return 1;
}

bool RenderCubes ()
{
    glutMainLoop ();
}
