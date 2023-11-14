# Verlet Cloth Mesh Component 

## A PBD Cloth Solver Plugin for Unreal Engine 4

<img src="Doc/gif_ClothVolumePreservation.gif" alt="Cloth Volume Preservation Demo GIF" style="zoom:50%;" />

#### NOTE: This project was developed in 2020-2021 for UE4. If I get time I do plan to refactor into UE5, this will require some changes, including utilising the RDG. 

This is a personal project from 2020-21 which contains an Unreal Engine 4 plugin that implements a velocity-less Verlet based Position Based Dynamics (PBD) cloth solver. The implementation is heavily based on the [Advanced Character Physics](https://www.cs.cmu.edu/afs/cs/academic/class/15462-s13/www/lec_slides/Jakobsen.pdf) Paper by Thomas Jakobsen [1]. This allows for a very high performance cloth solver, with some trade-offs such as timestep dependence. 

Cloth is discretized as a set of particle positions from some input cloth triangle mesh. Currently distance constraints are built based off the meshes edges and thus manifold and relatively good input topology is needed for clean results. Both Self and World collision constraints are implemented; self collisions utilise a spatial hash grid as described in [2] to only search nearby neighbour particles, while world collisions utilise the engines abstracted Physx or Chaos implementation to search for particle-world collisions by shape overlaps.

While distance constraints are constructed only once, self and world collision constraints are updated per tick. Both collision constraints utilise inequality constraints, which aim to be satisifed via projection in order to resolve themselves. Currently I implement a Gauss-Seidel like relaxation method for solving all constraints. Note that collision constraints are projected separately to the distance constraints, this is not accurate to the paper, but was done to reduce branching based on the currently enabled constraints. 
I most likely will remove this in-favour of a unified solver later, as it reduces convergence. 

This plugin utilises Epic's own Procedural Mesh Component with the 'VerletClothMeshComponent' derived from it. This allows us to efficiently update the cloth mesh on the GPU, without re-building per frame. In the future this could be used to allow tearing. All simulation is done within the game-thread, the current Gauss-Seidel relaxation method is inherently single threaded (unless a graph-colouring or Jacobi scheme is used instead). 

Particles do not explicitly store velocities, instead as shown in [1] the velocities are instead derived from current and previous positions $(\mathbf{x}_{n}, \mathbf{x}_{n-1})$. The velocity-less Verlet integrator directly computes $\mathbf{x}_{n+1}$ given the particles acceleration:
$$
\mathbf{x}_{n+1} = 2\mathbf{x}_n - \mathbf{x}_{n-1} + \ddot{\mathbf{x}}_n \Delta t^2
$$
This makes it efficient, while been second order accurate. In hindsight, I would not recommend this and instead use an alternative like Störmer-Verlet or Semi-Implicit Euler, as this allows control of velocities along with explicit storage. Furthermore, as all forces (beyond gravity) acting on the particles are typically resolved by constraints (within a PBD like framework), the integration accuracy is not super important. 

This is not the case here, as I implemented a volume-pressure conservation like force, which is resolved without constraints. Given a closed mesh, we can approximate its volume using average particle distances based on random sample positions, this will not work on an open mesh. 

Some key points: 

* This is a fast and simple solver, it could be great for set-dressing within the editor, by converting cloth back to static meshes. 
* I did write a GPU version utilising Compute shaders, I may refactor this as part of a port to UE5 and RDG. 
* This was a personal project from a few years ago, has not been maintained since me re-writing this README (due to other projects). 
* If I was writing a cloth solver plugin now I would not PBD + Verlet Integration, I would use XPBD, given it's stability and timestep invariance. 
* I've left this project online as it seems to have helped people writing plugins and understanding the aforementioned papers.  

## Usage

Simply place a empty actor into the scene, add a "Verlet Cloth Mesh" Component. Then select some Static Mesh to use as the base cloth mesh, and click the **BuildClothState** button, from here you can play with the simulation settings, while the solver ticks in the editor viewport, or run the component in PIE mode. 

## Building
The Plugin is was tested with  UE4.25 | 4.26 | 4.27

The plugin should be placed directly in a project location eg **./MyProject/Plugins**. This may require you to rebuild the module on your machine, this requires Visual Studio 2017/19/22 and the required components for UE4 C++ Development. 

The plugin comes with a few different test meshes, some of which are closed meshes. Note the Stanford Bunny asset is a re-meshed version of the Bunny from the Stanford 3D Scanning Library. 

## Known Issues
This project is a self educational project, and not intended for use in production, it has several known issues. 
* Self Collisions are resolved using discrete collision detection and are not completely stable. 
* If you change the Static Cloth Mesh within the component once **BuildClothState** has been invoked, the plugin will crash, this is due to the Static Mesh -> Procedural Mesh creation not been refreshed currently.

Nonetheless, this code is licensed under the MIT License, feel free to do what you wish with it. 

____

### References 
1. * T.Jakobsen, Advanced Character Physics, (IO Interactive).
2. * Particle-Based Fluid Simulation for Interactive Applications (Müller et al. SCA 2003).