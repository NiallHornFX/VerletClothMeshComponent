# Verlet Cloth Mesh Component Plugin - Implemented in Unreal Engine 4

![Cloth Patch Demo GIF](Doc/VerletClothPatch.gif)

## Implementation Overview 
The implementation is heavily based on the [Advanced Character Physics](https://www.cs.cmu.edu/afs/cs/academic/class/15462-s13/www/lec_slides/Jakobsen.pdf) Paper by Thomas Jakobsen [1] where cloth is discretized as a set of discrete particle positions from some input cloth manifold. In my case I chose to implement the ability to simulate cloth of arbitrary triangle mesh from a Static Mesh input, where the shared edges of vertices form the constraints connecting the cloth mesh together. These constraints are then projecting particle pair positions, to try and minimize the delta of there current distance back to the rest distance (derived from distance of the original vertices of the input static mesh) of the segment between them by some user stiffness coefficient k. It is also inspired by the Cable Component plugin which ships with Unreal Engine, that is solving a similar problem, but for a 1 Dimensional procedural cable. I Implemented this based on a subclassed Procedural Mesh Component, that both performs the simulation and rendering (via Procedural Mesh Sections) of the cloth mesh. 

The discretization is not unlike a simple Mass Spring based cloth model, however it differs in terms of how the equations of motion are solved, specially integration uses the velocity-less Verlet method oppose to something like the Forward (Explicit) Euler method. Unlike Euler where we use Force to integrate first derivative's (acceleration (a = f/m) to velocity, and velocity to position), Verlet can be seen as integrating second derivative's ie from Acceleration directly to position, thus velocity is implicitly calculated using the stored previous (x(n-1)) and current (x(n)) position to integrate the new (x(n+1)) position of each particle.

Also Oppose to using spring and damper forces to resolve constraints, we assume the particles are connected by springs of infinite stiffness and directly project there positions to satisfy constraints. However to satisfy constraint's of such a connected system, a Gauss-Seidel relaxation like approach is used, where constraints are projected iteratively, over multiple iterations to allow convergence to a (mostly) satisfied solution where the effects of constraints on particles on constraints (and so on) can be distributed throughout the system, minimizing the position delta's of pairs of particles along each constraint segment and thus creating cloth which can deform and crease as expected.  A similar approach is used for self collision's, however we are minimizing the intersection distance of particle radii instead of position deltas along triangle edge segments. For Self Collision I implemented a basic Spatial Hash Grid to accelerate the Particle-Particle distance and radii intersection checking code to remove the O(n^2) bottleneck of such an approach. The Spatial Hash function is based on the presentation from NVIDIA's Simon Green [NVIDIA's Simon Green talk](http://developer.download.nvidia.com/presentations/2008/GDC/GDC08_ParticleFluids.pdf).

I also added a bunch of Debug Functions to visualize the particle discretizations, the Hash Grid spatial locality on the particles and so forth inside the editor. 

## Building
This plugin is currently implemented as a project plugin, placed into a plugin directory eg **./MyProject/Plugins**. The plugin comes with a few different test meshes, some of which are closed meshes which
are for my future plans of implementing a fake volume preservation position projection solver step. Note the Stanford Bunny asset is from the Stanford 3D Scanning Library. 
The Plugin is currently built under Unreal Engine 4.25.1, however it should build with 4.26. 
Copy the contents into a folder named "VerletClothMesh" into your UE4 Project's Plugins directory, or your UE4 Engine Plugins install directory '../Epic Games/UE_4.2x/Engine/Plugins' and enable it in the editor. Note this may require you to rebuild the module on your machine, this requires Visual Studio 2017 and the required components for UE4 C++ Development. 
Then simply place a empty actor into the scene, add a "Verlet Cloth Mesh" Component, select some Static Mesh to use as the base cloth mesh, and click the "BuildClothState" button, from here you can play with
the simulation settings, while the solver ticks in the editor viewport, or run the component in PIE mode. 

## Known Issues
This project is a self educational project, I am still working on refinement's in a separate private repo -
* Performance with high resolution meshes is not ideal, theres no use of Multi-Threading or GPU Compute based acceleration yet. 
* The Tangents (and Normals) are not re-calculated per frame to update the procedural mesh sections yet. This will be fixed soon. 
* Self Collisions are still a work in progress, and no doubt will need to use a more complex projection step than currently implemented, to be more temporally stable.  
* Closed meshes are not correctly solved currently due to the lack of volumetric constraints, i'm working on a volume preservation method described briefly above. 
* If you change the Static Cloth Mesh within the component once BuildClothState has been invoked, the plugin will crash, this is due to the Static Mesh -> Procedural Mesh creation not been refreshed currently.

### References 
1. * T.Jakobsen, Advanced Character Physics, (IO Interactive).
2. * S.Green, Particle Based Fluid Simulation, (NVIDIA, GDC 2008).