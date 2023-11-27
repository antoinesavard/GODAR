# GODAR
## Granular Discrete Arctic Rheology: DEM for sea ice modeling

### What does GODAR do?

GODAR is a discrete element model that solves Newton's equations for a given number of cylindrical floes of ice. The disks have a radius, a thickness initial positions and bonds specified by the user, from which the model advances time. The forces at play are contact forces computed from Hertzian mechanics, bond forces computed from a simple beam theory, coriolis, atmopheric and oceanic forces computed from a quadratic drag law. Sheltering is taken into account, meaning that small floes can hide behind large floes and receive less wind or currents because of that. In order to find all interactions, a nearest neighbor search algorithm has been implemented: we use a kd-tree to find interacting floes. Moreover, the code is parallelized using mpi and openmp, and therefore GODAR is suitable for large number of floes on the order of $\mathcal{O}(10^3)$ particles.

### How to setup GODAR

There are a few things that needs to be done before you can compile and run this code. First off, the KdTree algorithm used in this program comes from coretran, so you need to install coretran on your machine. Coretran is available on Github at the following link: https://github.com/leonfoks/coretran. 

Once this is done, here are the steps you have to do. Let's call the path where you installed coretran: `/mycoretran`. There should be two subdirectories in your install folder: `/mycoretran/lib` and `/mycoretran/include`, in which you will find the `.so` file and the modules respectively.
We will need these in the SConstruct file; this is specific to your machine.

1. Change the `LIBPATH` and `F90PATH` values in the SConstruct file to your `/mycoretran/lib` and `/mycoretran/include` respectively
2. Next you will need to add the library path to your external environment variable. To do this we will edit the `~/.bashrc` file so that we only have to do this operation once, otherwise, you will have to do it everytime you open a new terminal. Open the terminal and run: `emacs ~/.bashrc`
3. Add the following line to this file: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mycoretran/lib` If the line already exist, just append `:/mycoretran/lib` at the end of the line.
4. `Crtl-x Crtl-s` to save and `Crtl-x Crtl-c` to exit emacs.
5. Then enter the `exec bash` to reload your terminal with the changes.

After these steps, you should be good to go. 

Next, to compile, this program uses `scons`. You will have to install it prior to running this program. The proper way of doing this is by creating a virtual environment first in which we will be able to install scons. If you already have your own python environment please skip ahead.

1. `python -m venv path/to/my/env` to create a virtual environment named `env`
2. The there are two options, you either always run the command `source path/to/my/env/bin/activate` or you put this line in your `~/.bashrc` file so that everytime you open a session, this specific environment gets loaded. To deactivate the environment: `deactivate`.
3. Then you install scons `python -m pip install scons`. 

After that, you need to execute the `init.sh` file in order to setup the folders. It will create an output folder and some default input files that you can edit. These files are actually in the `generic` folder, so that you can modify the ones created by the `init.sh` as you please (if you don't remember the syntax, you can go have a look at the ones in the `generic` directory).

Then, run `scons` with the appropriate parameters to compile Godar.

- To compile the code on n cores `scons -j n` 
- To clear the build: `scons -c`
- Debug the code: `scons debug=1`
- Run the executable in the background: `sh start.sh`

You can use the `start.sh` file if you are working on your own machin. It contains the following line:

```bash
mpirun --bind-to none -n 1 ./godar < input_restart > out
```

The `--bind-to none` makes all threads available to the program (you can then specify this number in the openmp part), while the `-n 1` simply tells how many copies of the program you want to run (on different machines for example). For example, let's say your machine has 2 cores 12 threads per core, you could either run 1 process on 24 threads using openmp or 2 processes (1 per core) with 12 threads using both mpi and openmp.
When working on clusters, you will usually encounter scheduler like slurm. In that case, don't use the `start.sh` file.

The input_file is a simple file to pass along to the main program when executing that can take the following arguments in this order: 
- NEED`[bool]   input namelist`
- NEED`[bool]   input restart`
- OPT `[int]    restart version`
- OPT `[int]    restart time`
- NEED`[int]    exp version`
- NEED`[int]    number of threads`

The code executes automatically a python script to output a video file of the collisions. This is saved in collisionXX.mp4 file. To suppress, just comment the line 

```fortran
call execute_command_line("python /plots/video.py")
```

To use it, change the path to your python interpreter and the location of the python code.

### How to modify the parameters in GODAR

All physical and numerical parameters as well as options are contained in the namelist.nml file. Therefore, if you want to play with the physics, you can simply change the values in this file without recompiling the code; this permits the use of batch job using an appropriate bash script that modifies the parameters that you want. If you forget the default parameters, they are all set in the `get_default` subroutine in the par_get.f90 file. You can also opt to not use the namelist by setting the read namelist option to false in the input file. 

The only parameters that are set at compilation are the number of particles (because we use this value for setting the lenght of all arrays), the size of the domain (for future developement where we will superimpose a grid over the domain to have spatially varying forcings), and the rank of the master thread. We are curently worknig on a way to work around this issue. (Probably using the input files and allocatable arrays in modules.)

The input files are setting the run informations: experiment number, whether to use the namelist or not, whether restarting from a previous experiement or not, etc. (This is why we think setting the particle number in there would be good.)

### How to improve GODAR

Godar is a very simple model. Here are some ideas to explore, or that we want to include in the future:

- spatially varying forcings (needs integrating spatially over the particle's area -> needs an underlying grid)
- 3D capabilities for waves (straightforward)
- fancier domain definition (needs a grid) for realistic Arctic domain
- rheological model inside the particle to define how particles can break into smaller particles (needs polygonal particles)
- polygonal particles (needs multicontact capabilities)

### Program architecture

```bash
- godar/
    - init.sh
    - generic/
        - input_genreric
        - input_restart_genreric
        - namelist_genreric.nml
        - SConstruct_genreric
        - start_genreric.sh
    - src/
        - bonds.f90
        - boundaries.f90
        - compaction.f90
        - contact_forces.f90
        - forcing.f90
        - ice.f90
        - ini_get.f90
        - integration.f90
        - mpi_coms.f90
        - par_get.f90
        - physics_bonds.f90
        - physics.f90
        - sea_ice_post.f90
        - stepper.f90
        - CB_bond.h
        - CB_const.h
        - CB_forcings.h
        - CB_variables.h
        - CB_mpi.h
        - CB_options
        - parameter.h
    - datetime/
        - date_mod.f90
        - datetime_mod.f90
        - datetimedelta_mod.f90
        - dateutils_mod.f90
        - precision_mod.f90
        - time_mod.f90
    - plots/
        - files.py
        - init.py
        - test.py
        - video.py
    - output/
    - input_files/
```