/*
Author: Alex Rodriguez
Class: ECE4122  
Last Date Modified: 09/19/2020
Description: This program develops a halftime show using UAVs for the final project.
Uses openGL and MPI to demo the show.
Program simulates UAVs flying up to a sphere and revolving around the sphere for ~60 seconds.
Used Red UAVs since my last name starts with R(our choice of color)
*/
#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "iomanip"
#include <cmath>
#include <math.h>
#include <cstdlib>
//#include <GLUT/glut.h>
#include <chrono>
#include <thread>
#include "ECE_Bitmap.h"

#define PI 3.14159
double deltaU = 1.0;
double u = 0.0166; // multiplied by deltaU for a timestep 
int onSurface = 0; // flag to check if on surface

const int numElements = 6; // for x, y, z, vx, vy, vz
const int rcvSize = 16 * 6; // (Main task + 15 UAVs) * numElements
//double* rcvbuffer = new double[rcvSize]; // receive buffer
double rcvbuffer[rcvSize]; // receive buffer on stack instead of heap in this case of known processes always being 16
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

// initializes positions for the first time
double initPosUAVs[45] = {0.0,-45.72,0.0, 0.0,45.72,0.0,23.46,-45.72,0.0,23.46,45.72,0.0,-23.46,-45.72,0.0,0.0,22.86,0.0,-23.46,22.86,0.0,23.46,22.86,0.0,
    -23.46,-22.86,0.0,23.46,-22.86,0.0,-23.46,0.0,0.0,23.46,0.0,0.0,0.0,0.0,0.0,0.0,-22.86,0.0,-23.46,45.72,0.0};

GLuint texture[1]; // will hold the texture for bitmap

BMP inBitmap; // bitmap


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

void processKey(unsigned char key, int x, int y)
{
	switch(key) 
	{
		case 'q':
		std::cout << "Done from rank: 0" << std::endl;
		//free(rcvbuffer);
		MPI_Abort(MPI_COMM_WORLD,0);
		//MPI_Finalize();
		//glutLeaveMainLoop();
		//exit(0);
		break;
	}
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
    glBindTexture(GL_TEXTURE_2D, 0);
}
/*
Function takes in no inputs and returns nothing
Displays each uav in their correct positions at the 0, 25, 50, 25, 0 yard-lines
Updates each position based on the receive buffer after the uavs take off
*/
void drawUAVs()
{
    //y coord: 54.864 - 18.288 = 36.576 -> 0 yard line
    //x coord: 23.468076 -> 0 yard line
    glColor3f(1.0, 0.0, 0.0); 
    for(int i = 0; i < 15; i++)
    {
    	glPushMatrix();
    		if(count == 0) // if uav hasnt taken off
    		{
    			glTranslatef(initPosUAVs[i*3], initPosUAVs[i*3+1], initPosUAVs[i*3+2]); // i =0 [0,1,2]// i =1 [3,4,5]//i=2 [6,7,8] // i=3 [9,10,11]// i=14 =[42,43,44]
    			glutSolidSphere(1, 20, 20); //1m sphere

    		}
    		else // draw in new position
    		{
    			glTranslatef(rcvbuffer[i*6+6], rcvbuffer[i*6+7], rcvbuffer[i*6+8]); //i=0 [6,7,8] // i =1 [12,13,14] i=2 [18,19,20]
            	glutSolidSphere(1, 20, 20); //1m sphere
    		}

    	glPopMatrix();

    }

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
        glutWireSphere(9.54, 10, 10); // wire sphere centered at 0,0,50 with 9.54m radius
    glPopMatrix();

    glutSwapBuffers(); // Make it all visible
    
    // gathers every rank info
    MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE, &rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);
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

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // Create Textures

    glGenTextures(1, texture);
    
    // Setup texture
    glBindTexture(GL_TEXTURE_2D, texture[0]);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); //scale linearly when image bigger than texture

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); //scale linearly when image smalled than texture


    glTexImage2D(GL_TEXTURE_2D, 0, 3, inBitmap.bmp_info_header.width, inBitmap.bmp_info_header.height, 0,
        GL_BGR_EXT, GL_UNSIGNED_BYTE, &inBitmap.data[0]);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL); // sets up environment of texture "color", no blending or anything like that

    glEnable(GL_TEXTURE_2D);

    glutDisplayFunc(renderScene); // renders all of scene
    glutReshapeFunc(changeSize); // changes size of window in case user expands window
    glutKeyboardFunc(processKey);
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
        uavForceX = 0.3; // 5N force
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
    else if(distP1P2 <= 12.0) // checks if close to sphere and ready to revolve around it
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
        double r = 1.0;//1.0;

        // parametric equations for calculating new positions every timestep around a sphere surface
        x = -r*sin(2*PI*(u*deltaU)); // changed from sin
        y = r*cos(2*PI*(u*deltaU)); // changed from cos
        z = 0; // changed from 0

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
        sendBuffer[0] = x+ rcvbuffer[rank*5+rank]; // new x position
        sendBuffer[1] = y+ rcvbuffer[rank*5+rank+1]; // new y position
        sendBuffer[2] = z+ rcvbuffer[rank*5+rank+2]; // new z position
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
        std::cout << "got out of mainloop!" << std::endl;
    	//MPI_Finalize();
    }
    else
    {
        // Sleep for 5 seconds
        std::this_thread::sleep_for(std::chrono::seconds(5));
        for (int ii = 0; ii < 600 ; ii++) // runs for about 68 seconds
        {
        	// Sleep for 5 seconds
            ++flag;
            CalcualteUAVsLocation(rank); // each rank handles a uav
            MPI_Allgather(sendBuffer, numElements, MPI_DOUBLE, &rcvbuffer, numElements, MPI_DOUBLE, MPI_COMM_WORLD);
        }
        std::cout << "Done from rank:  " << rank << std::endl;
        {
        	// free(rcvbuffer);
        	// rcvbuffer = NULL;
        	// std::cout << "Freed from rank: " << rank << std::endl;
        	MPI_Abort(MPI_COMM_WORLD,1);

        }
    }

    return 0;
}