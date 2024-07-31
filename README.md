# GODAR

## Granular Discrete Arctic Rheology: DEM for sea ice modeling

### What does GODAR do?

GODAR is a discrete element model that solves Newton's equations for a given number of cylindrical floes of ice. The disks have a radius, a thickness initial positions and bonds specified by the user, from which the model advances time. The forces at play are contact forces computed from Hertzian mechanics, bond forces computed from a simple beam theory, coriolis, atmopheric and oceanic forces computed from a quadratic drag law. Sheltering is taken into account, meaning that small floes can hide behind large floes and receive less wind or currents because of that. In order to find all interactions, a nearest neighbor search algorithm has been implemented: we use a kd-tree to find interacting floes. Moreover, the code is parallelized using mpi and openmp, and therefore GODAR is suitable for large number of floes on the order of $\mathcal{O}(10^3)$ particles.

### How to setup GODAR

#### Creating the config files

The first step will be to execute `./init.sh` file in order to setup the folders. It will create an output folder and some default input files that you can edit. These files are actually in the `generic` folder, so that you can modify the ones created by the `init.sh` as you please (if you don't remember the syntax, you can go have a look at the ones in the `generic` directory).

#### Intalling coretran

Next you will need to install a few things. There are a few things that needs to be done before you can compile and run this code. First off, the KdTree algorithm used in this program comes from coretran, so you need to install coretran on your machine. Coretran is available on Github at the following link: <https://github.com/leonfoks/coretran>.

I would suggest to follow the detailed instructions provided in coretran's readme as it well written and easy to use. I would recommend that coretran be installed in:

OSX: /opt/coretran

Linux: /usr/local/lib

Once this is done, here are the steps you have to do. Let's call the path where you installed coretran: `/mycoretran`. There should be two subdirectories in your install folder: `/mycoretran/lib` and `/mycoretran/include`, in which you will find the `.so` file (`.dylib` on OSX) and the modules respectively.
We will need these in the SConstruct file; this is specific to your machine.

##### Linux

The following fix should work on all machines, although it is a bit of a work-around.

1. Change the `LIBPATH` and `F90PATH` values in the SConstruct file to your `/mycoretran/lib` and `/mycoretran/include` respectively
2. Next you will need to add the library path to your external environment variable. To do this we will edit the `~/.bashrc` file so that we only have to do this operation once, otherwise, you will have to do it everytime you open a new terminal. Open the terminal and run: `emacs ~/.bashrc`
3. Add the following line to this file: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mycoretran/lib` If the line already exist, just append `:/mycoretran/lib` at the end of the line.
4. `Crtl-x Crtl-s` to save and `Crtl-x Crtl-c` to exit emacs.
5. Then enter the `exec bash` to reload your terminal with the changes.

##### OSX

On mac, there is an extra step to be taken as well.

1. Change the `LIBPATH` and `F90PATH` values in the SConstruct file to your `/mycoretran/lib` and `/mycoretran/include` respectively
2. In the `start.sh` file, copy the following line above the `mpirun` command. This line will tell the executable where to look for the dynamically shared libraries.

```bash
install_name_tool -change @rpath/libcoretran.dylib /mycoretran/libcoretran.dylib godar
```

You will also have to install a few things using brew. First, install brew as said on their website <https://brew.sh>.

Then you will need to install `gcc`, `open-mpi`, and `libomp`.

When using `mpifort` with brew, there might be an issue in compatibilities where you get an error message saying that the linker `ld` is not able to find the flag `-ld_classic`. If that is the case, then you have to move to the folder where your `homebrew` is and change the following file:

```bash
/opt/homebrew/share/openmpi/mpifort-wrapper-data.txt
```

In this file, there is a line that starts with `compiler_flags=` and you have to remove everything except the first part, so that it looks something like:

```bash
compiler_flags=-I${includedir}
```

This line is essentially just telling which flags to use when compiling the code with mpifort. The thing is that with xcode 15+, these flags are now obsolete, and are not compatible anymore with brew, well, as least by my findings. Removing these flags fixed things for me.

After these steps, you should be good to go.

#### Installing scons and creating a virtual environment

To install packages on a super computer, it is usually better to do that in a virtual environment. If you need a virtual environment, create one first. If not, skip to the installation of scons.

##### Optional virtual environment

1. `python -m venv path/to/my/env` to create a virtual environment named `env`
2. The there are two options, you either always run the command `source path/to/my/env/bin/activate` or you put this line in your `~/.bashrc` file so that everytime you open a session, this specific environment gets loaded. To deactivate the environment: `deactivate`.

##### Install scons

Next, to compile this program, we use `scons`. You will have to install it prior to running this program. The proper way of doing this is by creating a virtual environment first in which we will be able to install scons. If you already have your own python environment, or do not need one, then you can proceed.

You install scons via pip or conda if that's what you use: `python -m pip install scons`. Here are some useful commands of scons.

- To compile the code on n cores `scons -j n`
- To clear the build: `scons -c`
- Debug the code: `scons debug=1`
- Run the executable in the background from GODAR/jobs/: `./start.sh`

#### How to run godar

After compilation, you can use the `start.sh` file if you are working on your own machine. It contains the following line:

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

### How to modify the parameters in GODAR

All physical and numerical parameters as well as options are contained in the namelist.nml file. Therefore, if you want to play with the physics, you can simply change the values in this file without recompiling the code; this permits the use of batch job using an appropriate bash script that modifies the parameters that you want. If you forget the default parameters, they are all set in the `get_default` subroutine in the par_get.f90 file. You can also opt to not use the namelist by setting the read namelist option to false in the input file.

The only parameters that are set at compilation are the number of particles (because we use this value for setting the lenght of all arrays), the size of the domain (for future developement where we will superimpose a grid over the domain to have spatially varying forcings), and the rank of the master thread. We are curently worknig on a way to work around this issue. (Probably using the input files and allocatable arrays in modules.)

The input files are setting the run informations: experiment number, whether to use the namelist or not, whether restarting from a previous experiement or not, etc. (This is why we think setting the particle number in there would be good.)

### How to run large number of simulations using SLURM

Usually, you will run GODAR on a super computer, which requires a job scheduler like slurm. The scripts in /utils have been written with that usage in mind.

First, you can run `./run_init.py init_args.dat`, which will create initial conditions for your runs. It reads the input_init.dat file, where you can enter various distributions and parameters for both radius and thickness, the rest is handled by the program.

Second, you can run `./duplication.sh < args.dat` which will create namelists based on the arguments provided in the input file. This is for when you want to run the same experiment with multiple types forcings, or different experiments with the same forcing.

Third, you simply run the slurm file the previous program created in `/jobs`. This is where you can modify the allocated time and memory, etc.

Finally, when your jobs have finished running, you can run `./run_video.sh video_args.dat` to create all the videos.

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
    - datetime/
        - date_mod.f90
        - datetime_mod.f90
        - datetimedelta_mod.f90
        - dateutils_mod.f90
        - precision_mod.f90
        - time_mod.f90
    - files/
    - generic/
        - input_genreric
        - input_restart_genreric
        - namelist_genreric.nml
        - SConstruct_genreric
        - start_genreric.sh
    - inputs/
    - jobs/
    - namelist/
    - output/
    - plots/
        - anim/
        - plot/
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
        - physics.f90
        - reset.f90
        - sea_ice_post.f90
        - stepper.f90
        - CB_bond.h
        - CB_const.h
        - CB_forcings.h
        - CB_variables.h
        - CB_mpi.h
        - CB_options
        - parameter.h
    - tools/
        - analysis/
            - monte_carlo.py
        - init/
            - duplication_2nd_plate.py
            - gridded_placement.py
            - init2.py
        - plotter/
            - figures.py
            - video.py
            - void_ratios.py
        - tests/
            - scaling.py
            - test.py
        - utils/
            - files.py
    - utils/
        - duplication.sh
        - run_init.sh
        - run_video.sh
```
