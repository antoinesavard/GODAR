# GODAR
## Granular Discrete Arctic Rheology: DEM for sea ice modeling

### What does GODAR do?

GODAR is a discrete element model that solves Newton's equations for a given number of cylindrical floes of ice. The disks have a radius, a thickness initial positions and bonds specified by the user, from which the model advances time. The forces at play are contact forces computed from Hertzian mechanics, bond forces computed from a simple beam theory, coriolis, atmopheric and oceanic forces computed from a quadratic drag law. Sheltering is taken into account, meaning that small floes can hide behind large floes and receive less wind or currents because of that. In order to find all interactions, a nearest neighbor search algorithm has been implemented: we use a kd-tree to find interacting floes. Moreover, the code is parallelized, and therfore GODAR is suitable for large number of floes $\mathcal{O}(10^3)$.

### How to setup GODAR

There are a few things that needs to be done before you can compile and run this code. First off, the KdTree algorithm used in this program comes from coretran, so you need to install coretran on your machine. Coretran is available on Github. Once this is done, here are the steps you have to do. Let's call the path where you installed coretran: `/mycoretran`. There should be two subdirectories in your install folder: `/mycoretran/lib` and `/mycoretran/include`.

1. Change the `LIBPATH` and `F90PATH` values in the SConstruct file (lines 62 and 63) to your `/mycoretran/lib` and `/mycoretran/include` respectively
2. Open the terminal and run: `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mycoretran/lib`
3. Then run the command: `export LD_LIBRARY_PATH`

After these steps, you should be good to go. Next, to compile, run `scons-3` with the appropriate parameters.

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
call execute_command_line("python /plots/video.py")
```

To use it, change the path to your python interpreter and the location of the python code.
