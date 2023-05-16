# GODAR
## Granular Discrete Arctic Rheology: DEM for sea ice modeling

DEM model. To compile, run scons-3 with the appropriate parameters.

- To compile the code on n cores `scons-3 -j n` 
- To clear the build: `scons-3 -c`
- Debug the code: `scons-3 debug=1`
- Run the executable: `./godar < input_file`

The input_file is a simple file to pass along to the main program when executing that can take the following arguments in this order: 
- NEED`[bool]   input namelist`
- NEED`[bool]   input restart`
- OPT `[int]    restart version`
- OPT `[int]    restart time`
- NEED`[int]    exp version`
- NEED`[int]    number of threads`

The code executes automatically a python script to output a video file of the collisions. This is saved in collisionXY.mp4 file. To suppress, just comment the line 

```fortran
call execute_command_line("/aos/home/asavard/anaconda3/bin/python /storage/asavard/DEM/plots/video.py")
```

To use it, change the path to your python interpreter and the location of the python code.

What does GODAR do?
---
GODAR is a discrete element model that solves Newton's equations for a given number of cylindrical floes of ice. The disks have a radius, a thickness initial positions and bonds specified by the user, from which the model advances time. The forces at play are contact forces computed from Hertzian mechanics, bond forces computed from a simple beam theory, coriolis, atmopheric and oceanic forces computed from a quadratic drag law. Sheltering is taken into account, meaning that small floes can hide behind large floes and receive less wind or currents because of that. In order to find all interactions, a nearest neighbor search algorithm has been implemented: we use a kd-tree to find interacting floes. Moreover, the code is parallelized, and therfore GODAR is suitable for large number of floes $\mathcal{O}(10^3)$.
