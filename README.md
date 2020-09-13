# 3D-Half-Time-Show-Simulation
 A 3D simulation using MPI and OpenGL to render a half-time show on a football field with UAVs

Implements a 3D rendering of 15 UAVs (drawn as red glut sphere shapes) on a football field, with 1kg masses, a 1m cube bounding box, and able to generate a single force vector with a total magnitude of 20N in any direction. Simulates UAV physics with gravity when UAVs fly toward a virtual sphere in the sky. Once UAVs have reached the surface of the sphere, i.e 10 meters to the center point of the sphere, they fly along the surface for 60 seconds.

The simulation is implemented using an MPI application using 16 processes. The main rank is responsible for rendering the 3D scene with the 15 UAVs and football field while the other 15 processes are each responsible for controlling the motion of a single UAV.

The main process gathers the location and velocity vector of each UAV every 100ms and broadcasts this information to all the UAV processes and updates the 3D scene.
The flight of the UAV is controlled by normal kinematic rules.
