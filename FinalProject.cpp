/*
Author: Alex Rodriguez
Class: ECE4122  
Last Date Modified: 12/3/2019
Description: This program develops a halftime show using UAVs for the final project.
Uses openGL and MPI to demo the show.
Program simulates UAVs flying up to a sphere and revolving around the sphere for ~60 seconds.
Used Red UAVs since my last name starts with R(our choice of color)
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "iomanip"
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <GL/glut.h>
#include <chrono>
#include <thread>
#include "ECE_Bitmap.h"

#define PI 3.14159
double deltaU = 1.0;
double u = 0.0166; // multiplied by deltaU for a timestep 
int onSurface = 0; // flag to check if on surface

const int numElements = 6; // for x, y, z, vx, vy, vz
const int rcvSize = 16 * 6; // (Main task + 15 UAVs) * numElements
double* rcvbuffer = new double[rcvSize]; // receive buffer
double sendBuffer[numElements]; // send buffer

int volatile flag = 0; // check if positions have been drawn at first
int count = 0;
double distP1P2; // distance from one point to another

double velMag = 0; // magnitude of velocity aka speed
double velX = 0; // velocity in x direction
double velY = 0; // velocity in y direction
double velZ = 0; // velocity in z direction
double mass = 1; // mass of a uav
double uavForceX; // force x for a UAV
double fX; // force uav x direction
double fY; // force uav y direction
double uavForceZ; // force z for a UAV
double fZ; // force uav z direction with gravity included

double initPosUAVs[45]; // intial positions of each UAV

GLuint texture[2]; // will hold the texture for bitmap

// texture image for size and data
struct Image {

    unsigned long sizeX;

    unsigned long sizeY;

    char *data;

};

typedef struct Image Image;

BMP inBitmap; // bitmap


#define checkImageWidth 64

#define checkImageHeight 64


GLubyte checkImage[checkImageWidth][checkImageHeight][3];

// checks image for texture
void makeCheckImage(void) {

    int i, j, c;

    for (i = 0; i < checkImageWidth; i++) {

        for (j = 0; j < checkImageHeight; j++) {

            c = ((((i & 0x8) == 0) ^ ((j & 0x8) == 0))) * 255;

            checkImage[i][j][0] = (GLubyte)c;

            checkImage[i][j][1] = (GLubyte)c;

            checkImage[i][j][2] = (GLubyte)c;

        }

    }

}


//----------------------------------------------------------------------
// Reshape callback
// Function takes in a width w and height h. Returns nothing
// Window size has been set/changed to w by h pixels. Set the camera
// perspective to 45 degree vertical field of view, a window aspect
// ratio of w/h, a near clipping plane at depth 1, and a far clipping
// plane at depth 100. The viewport is the entire window.
//
//----------------------------------------------------------------------
void changeSize(int w, int h)
{
    float ratio = ((float)w) / ((float)h); // window aspect ratio
    glMatrixMode(GL_PROJECTION); // projection matrix is active
    glLoadIdentity(); // reset the projection
    gluPerspective(60.0, ratio, 0.1, 1000.0); // perspective transformation
    glMatrixMode(GL_MODELVIEW); // return to modelview mode
    glViewport(0, 0, w, h); // set viewport (drawing area) to entire window
}


/*
Function takes in no input and returns nothing
Displays the football field texture on a rectangle with correct meter lengths and widths, centered at 0,0,0
*/
void displayFootballField()
{
    // field is 120 yards long w/ end zones, 109.728m -> /2 = 54.864
    // field is 100 yards long w/o end zones, 91.44m -> /2 = 45.72
    // field is 53 1/3 yards wide, 46.936152m -> /2 = 23.468076
    // x axis is width of field, y axis is length of field
    //3.0 added to adjust for bmp bounds in image
    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, texture[0]);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);
    glVertex2f(-(23.468+3.0), (54.864+3.0));
    glTexCoord2f(0, 1);
    glVertex2f((23.468+3.0), (54.864+3.0));
    glTexCoord2f(1, 1);
    glVertex2f((23.468+3.0), -(54.864+3.0));
    glTexCoord2f(1, 0);
    glVertex2f(-(23.468+3.0), -(54.864+3.0));
    glEnd();
    glPopMatrix();

}
/*
Function takes in no inputs and returns nothing
Displays each uav in their correct positions at the 0, 25, 50, 25, 0 yard-lines
Updates each position based on the receive buffer after the uavs take off
*/
void drawUAVs()
{
    glBindTexture(GL_TEXTURE_2D,0);
    //y coord: 54.864 - 18.288 = 36.576 -> 0 yard line
    //x coord: 23.468076 -> 0 yard line
    glColor3f(1.0, 0.0, 0.0); 
    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[0], initPosUAVs[1], initPosUAVs[2]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[6], rcvbuffer[7], rcvbuffer[8]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[3], initPosUAVs[4], initPosUAVs[5]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[12], rcvbuffer[13], rcvbuffer[14]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[6], initPosUAVs[7], initPosUAVs[8]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[18], rcvbuffer[19], rcvbuffer[20]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[9], initPosUAVs[10], initPosUAVs[11]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[24], rcvbuffer[25], rcvbuffer[26]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[12], initPosUAVs[13], initPosUAVs[14]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[30], rcvbuffer[31], rcvbuffer[32]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[15], initPosUAVs[16], initPosUAVs[17]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[36], rcvbuffer[37], rcvbuffer[38]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[18], initPosUAVs[19], initPosUAVs[20]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[42], rcvbuffer[43], rcvbuffer[44]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[21], initPosUAVs[22], initPosUAVs[23]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[48], rcvbuffer[49], rcvbuffer[50]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[24], initPosUAVs[25], initPosUAVs[26]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[54], rcvbuffer[55], rcvbuffer[56]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[27], initPosUAVs[28], initPosUAVs[29]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[60], rcvbuffer[61], rcvbuffer[62]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[30], initPosUAVs[31], initPosUAVs[32]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[66], rcvbuffer[67], rcvbuffer[68]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[33], initPosUAVs[34], initPosUAVs[35]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[72], rcvbuffer[73], rcvbuffer[74]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[36], initPosUAVs[37], initPosUAVs[38]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[78], rcvbuffer[79], rcvbuffer[80]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[39], initPosUAVs[40], initPosUAVs[41]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[84], rcvbuffer[85], rcvbuffer[86]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    glPushMatrix();
        if(count == 0) // if uav hasnt taken off
        {
            glTranslatef(initPosUAVs[42], initPosUAVs[43], initPosUAVs[44]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
        else // draw in new position
        {
            glTranslatef(rcvbuffer[90], rcvbuffer[91], rcvbuffer[92]);
            glutSolidSphere(1, 20, 20); //1m sphere
        }
    glPopMatrix();

    count++; // done with initial positions
}

//----------------------------------------------------------------------
// Draws the entire scene
// Function takes in no inputs and returns nothing
// We first update the camera location based on its distance from the
// origin and its direction.
// Calls all drawing functions from above
//----------------------------------------------------------------------
void renderScene()
{

    // Clear color and depth buffers
    glClearColor(0.0, 0.0, 0.0, 1.0); // background color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Reset transformations
    glLoadIdentity();

    // positions camera to see entirety of scene
    gluLookAt(90, 0, 100, 
              0, 0, 50,
              0.0, 0.0, 1.0);

    // loads correct model view
    glMatrixMode(GL_MODELVIEW);

    glPushMatrix();
    displayFootballField(); // calls football display function
    glPopMatrix();

    glPushMatrix();
    drawUAVs(); // calls uav display function
    glPopMatrix();

    glPushMatrix();
        glColor3f(0.0, 0.0, 1.0);
        glTranslatef(0, 0, 50.0);
        glutWireSphere(10, 10, 10); // wire sphere centered at 0,0,50 with 10m radius
    glPopMatrix();

    glutSwapBuffers(); // Make it all visible
    
    // gathers every rank info
    MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE, rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);
}

//----------------------------------------------------------------------
// timerFunction  - called whenever the timer fires
// Function takes in an int id and returns nothing
// calls itself every 100ms for a screen refresh rate 
//----------------------------------------------------------------------
void timerFunction(int id)
{
    glutPostRedisplay(); // redisplays scene to update anything that was changed every 100ms
    glutTimerFunc(100, timerFunction, 0);
}

//----------------------------------------------------------------------
// mainOpenGL  - standard GLUT initializations and callbacks
// Function takes in two ints and returns nothing
// calls all callbacks and sets up texture bitmap
//----------------------------------------------------------------------
void mainOpenGL(int argc, char**argv)
{
    glutInit(&argc, argv); // initializes gl
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100); 
    glutInitWindowSize(400, 400); // window size 400x400

    glutCreateWindow(argv[0]);
    glEnable(GL_DEPTH_TEST); // enables depth
    glClearColor(0.0, 0.0, 0.0, 0.0); // clears color
    glShadeModel(GL_SMOOTH);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    inBitmap.read("ff.bmp"); // reads in bitmap of football field

    makeCheckImage();

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // Create Textures

    glGenTextures(2, texture);
    
    // Setup first texture
    glBindTexture(GL_TEXTURE_2D, texture[0]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); //scale linearly when image bigger than texture

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); //scale linearly when image smalled than texture


    glTexImage2D(GL_TEXTURE_2D, 0, 3, inBitmap.bmp_info_header.width, inBitmap.bmp_info_header.height, 0,
        GL_BGR_EXT, GL_UNSIGNED_BYTE, &inBitmap.data[0]);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

    //Do the second texture

    glBindTexture(GL_TEXTURE_2D, texture[1]);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

    glTexImage2D(GL_TEXTURE_2D, 0, 3, checkImageWidth, checkImageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, &checkImage[0][0][0]);

    glEnable(GL_TEXTURE_2D);


    glutReshapeFunc(changeSize); // changes size of window in case user expands window
    glutDisplayFunc(renderScene); // renders all of scene
    glutTimerFunc(100, timerFunction, 0); // 100ms refresh rate called
    glutMainLoop();
}

/*
Function takes in a rank and returns nothing
Using the rank, controls each uav assigned to each rank
Updates position of uav using kinematics
Moves to center of sphere and then revolves around the surface of the sphere
*/

void CalcualteUAVsLocation(int rank)
{
    if(flag == 1)
    {
        //initialize stuff at first to be zero
        rcvbuffer[rank*5+rank+3] = 0; //vx0
        rcvbuffer[rank*5+rank+4] = 0; //vy0
        rcvbuffer[rank*5+rank+5] = 0; //vz0

        rcvbuffer[rank*5+rank] = initPosUAVs[rank*3-3]; // init x
        rcvbuffer[rank*5+rank+1] = initPosUAVs[rank*3-2]; // init y
        rcvbuffer[rank*5+rank+2] = initPosUAVs[rank*3-1]; // init z
    }
    //distance formula from uav to center of sphere
    distP1P2 = sqrt(pow((rcvbuffer[rank*5+rank] - 0), 2.0)+pow((rcvbuffer[rank*5+rank+1] - 0), 2.0)+pow((rcvbuffer[rank*5+rank+2] - 50), 2.0));
    // if distance from uav to point is greater than 11, fly towards it.
    if(distP1P2 > 11.0)
    {

        //generate a single force vector with a total magnitude of 20 N in any direction.
        // use up to 20 n of force, mainly in the z direction
        // direction vector from uav to center of sphere
        double dirX = 0.0 - rcvbuffer[rank*5+rank];
        double dirY = 0.0 - rcvbuffer[rank*5+rank+1];
        double dirZ = 50.0 - rcvbuffer[rank*5+rank+2];
        double magnitudeXYZ = sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ); // magnitude of direction vector
        // normalized direction vector
        double normDirX = dirX/magnitudeXYZ;
        double normDirY = dirY/magnitudeXYZ;
        double normDirZ = dirZ/magnitudeXYZ;
        double uavForceX = 0.5; // 5N force
        fX = 0.5; //(uavForceX*rank); ((double)rank)/3.0 // force in x direction
        fY = (8.0-(double)rank/2); // force in Y direction
        uavForceZ = 9.4; // base force needed to counter gravity
        fZ = (1+uavForceZ+(rank*0.5)-9.8); // Z force to counter act force of gravity, ranges from .10 to 7.1 after gravity is included
        double fMag = sqrt(fX*fX + fY*fY + fZ*fZ); // always less than 20N-9.8N

        if(fMag > 20.0) // if mag exceeds 20N, cut down force
        {
            fX = 2.0;
            fY = 2.0;
            fZ = 10.0;
        }

        double accelX = fX/mass; // acceleration in x direction
        double accelY = fY/mass; // acceleration in y direction
        double accelZ = fZ/mass; // acceleration in z direction

        // kinematic equations for x, y, z
        double newX = rcvbuffer[rank*5+rank] + rcvbuffer[rank*5+rank+3]*normDirX*0.1 + (0.5)*accelX*normDirX*(0.1*0.1);
        double newY = rcvbuffer[rank*5+rank+1] + rcvbuffer[rank*5+rank+4]*normDirY*0.1 + (0.5)*accelY*normDirY*(0.1*0.1);
        double newZ = rcvbuffer[rank*5+rank+2] + rcvbuffer[rank*5+rank+5]*normDirZ*0.1 + (0.5)*accelZ*normDirZ*(0.1*0.1);

        //kinematic equations for velocity in x,y,z direction
        double newVelX = rcvbuffer[rank*5+rank+3] + accelX*(0.1);
        double newVelY = rcvbuffer[rank*5+rank+4] + accelY*(0.1);
        double newVelZ = rcvbuffer[rank*5+rank+5] + accelZ*(0.1);

        double newVelMag = sqrt(newVelX*newVelX + newVelY*newVelY + newVelZ*newVelZ); // velocity magnitude
        if(newVelMag > 10.0) // cannot exceed 10m/s
        {
            // if exceeded 10m/s, apply force in negative direction to slow down
            fX = -accelX*mass;
            fY = -accelY*mass;
            fZ = -accelZ*mass;
        }


        sendBuffer[0] = newX; // new x position
        sendBuffer[1] = newY; // new y position
        sendBuffer[2] = newZ; // new z position

        sendBuffer[3] = newVelX; // new vx position
        sendBuffer[4] = newVelY; // new vy position
        sendBuffer[5] = newVelZ; // new vz position
    }
    else if(distP1P2 <= 11.0) // checks if close to sphere and ready to revolve around it
    {
        onSurface = 1; // is now on surface
        if(deltaU > 61) // timestep reset
        {
            deltaU = 1.0;
        }
        // initialized positions and radius
        double x = 0; 
        double y = 0;
        double z = 0;
        double r = 1;

        // parametric equations for calculating new positions every timestep around a sphere surface
        x = -r*sin(2*PI*(u*deltaU));
        y = r*cos(2*PI*(u*deltaU));
        z = 0;

        // new direction vector
        double dirX = 0 - x;
        double dirY = 0 - y;
        double dirZ = 50 - rcvbuffer[rank*5+rank+2];
        double magnitudeXYZ = sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ); // magnitude of direction vector
        // normalized direction vector
        double normDirX = dirX/magnitudeXYZ;
        double normDirY = dirY/magnitudeXYZ;
        double normDirZ = dirZ/magnitudeXYZ;

        // uses theta to differentiate each UAV orbit and recalculates x,y,z based on it
        double theta = (((double)rank-1.0)*12)*(PI/180);
        double xPrime = x;
        double yPrime = y*cos(theta) - z*sin(theta);
        double zPrime = y*sin(theta) + z*cos(theta);

        x = xPrime;
        y = yPrime;
        z = zPrime;
        ++deltaU; // increased time step for next iteration
        sendBuffer[0] = x + rcvbuffer[rank*5+rank]; // new x position
        sendBuffer[1] = y + rcvbuffer[rank*5+rank+1]; // new y position
        sendBuffer[2] = z + rcvbuffer[rank*5+rank+2]; // new z position
    }

}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Main entry point determines rank of the process and follows the 
// correct program path
// Function takes in an int and char and returns an int
// Function main handles MPI and gives rank 0 the rendering job and 
// rank 1-15 the update uav location job for each uav 1-15
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main(int argc, char**argv)

{
    // initializes positions for the first time
    initPosUAVs[0] = 0.0;
    initPosUAVs[1] = -45.72;
    initPosUAVs[2] = 0.0;
    initPosUAVs[3] = 0.0;
    initPosUAVs[4] = 45.72;
    initPosUAVs[5] = 0.0;
    initPosUAVs[6] = 23.46;
    initPosUAVs[7] = -45.72;
    initPosUAVs[8] = 0.0;
    initPosUAVs[9] = 23.46;
    initPosUAVs[10] = 45.72;
    initPosUAVs[11] = 0.0;
    initPosUAVs[12] = -23.46;
    initPosUAVs[13] = -45.72;
    initPosUAVs[14] = 0.0;
    initPosUAVs[15] = 0.0;
    initPosUAVs[16] = 22.86;
    initPosUAVs[17] = 0.0;
    initPosUAVs[18] = -23.46;
    initPosUAVs[19] = 22.86;
    initPosUAVs[20] = 0.0;
    initPosUAVs[21] = 23.46;
    initPosUAVs[22] = 22.86;
    initPosUAVs[23] = 0.0;
    initPosUAVs[24] = -23.46;
    initPosUAVs[25] = -22.86;
    initPosUAVs[26] = 0.0;
    initPosUAVs[27] = 23.46; //
    initPosUAVs[28] = -22.86;
    initPosUAVs[29] = 0.0;
    initPosUAVs[30] = -23.46;
    initPosUAVs[31] = 0.0;
    initPosUAVs[32] = 0.0;
    initPosUAVs[33] = 23.46;
    initPosUAVs[34] = 0.0;
    initPosUAVs[35] = 0.0;
    initPosUAVs[36] = 0.0;
    initPosUAVs[37] = 0.0;
    initPosUAVs[38] = 0.0;
    initPosUAVs[39] = 0.0;
    initPosUAVs[40] = -22.86;
    initPosUAVs[41] = 0.0;
    initPosUAVs[42] = -23.46;
    initPosUAVs[43] = 45.72;
    initPosUAVs[44] = 0.0;


    int numTasks;
    int rank;

    int rc = MPI_Init(&argc, &argv);

    if (rc != MPI_SUCCESS) 
    {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int gsize = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &gsize);

    if (rank == 0) 
    {
        mainOpenGL(argc, argv);
    }
    else
    {
        // Sleep for 5 seconds
        std::this_thread::sleep_for(std::chrono::seconds(5));
        for (int ii = 0; ii < 600 ; ii++) // runs for about 68 seconds
        {
            ++flag;
            CalcualteUAVsLocation(rank); // each rank handles a uav
            MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE, rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);
        }
    }

    //MPI_Finalize();
    return 0;
}