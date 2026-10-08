# GODAR

## Granular flOe Dynamics for seA ice Rheology

### What does GODAR do?

GODAR is a discrete element model that solves Newton's equations for a given number of cylindrical floes of ice. The disks have a radius, a thickness initial positions and bonds specified by the user, from which the model advances time. The forces at play are contact forces computed from Hertzian mechanics, bond forces computed from a simple beam theory, coriolis, atmopheric and oceanic forces computed from a quadratic drag law. A simple ridging scheme is provided. Sheltering is taken into account, meaning that small floes can hide behind large floes and receive less wind or currents because of that. In order to find all interactions, a nearest neighbor search algorithm has been implemented: we use a kd-tree to find interacting floes. Moreover, the code is parallelized using mpi and openmp, and therefore GODAR is suitable for large number of floes on the order of $\mathcal{O}(10^3)$ particles.

### How to setup GODAR

#### Installing basic dependencies

You will need cmake to compile the different libraries and source files.

##### OSX

You will also have to install a few things using brew. First, install brew as said on their website <https://brew.sh>.

Then you will need to install `gcc`, `open-mpi`, `ffmpeg`, and `libomp`.

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

#### Intalling coretran

Next you will need to install a few things. There are a few things that needs to be done before you can compile and run this code. First off, the KdTree algorithm used in this program comes from coretran, so you need to install coretran on your machine. Coretran is available on Github at the following link: <https://github.com/leonfoks/coretran>.

I would suggest to follow the detailed instructions provided in coretran's readme as it is well written and easy to use. I would recommend that coretran be installed close to where `GODAR/` is, e.g.: `GODAR/../coretran`. To do this, enter the following commands:

```bash
git clone git@github.com:leonfoks/coretran.git
cd coretran
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/GODAR/../coretran ../src
make install
```

Once this is done, here are the steps you have to do. Let's call the path where you installed coretran: `/mycoretran`. There should be two subdirectories in your install folder: `/mycoretran/lib` and `/mycoretran/include`, in which you will find the `.so` file (`.dylib` on OSX) and the modules respectively.
We will need these in the CMakeLists.txt file; this is specific to your machine, specifically the `/mycoretran/lib/cmake` directory. Cmake and the `CMakeLists.txt` file should deal with this automatically, but it is good for you to know.

#### Creating a virtual environment

To install packages on a super computer, it is usually better to do that in a virtual environment. If you need a virtual environment, create one first.

##### Optional virtual environment

1. `python -m venv path/to/my/env` to create a virtual environment named `env`
2. The there are two options, you either always run the command `source path/to/my/env/bin/activate` or you put this line in your `~/.bashrc` file so that everytime you open a session, this specific environment gets loaded. To deactivate the environment: `deactivate`.

#### Creating the config files and installing GODAR

The first step will be to execute `./utils/init.sh` file in order to setup the folders. It will create an output folder and some default input files that you can edit. These files are actually in the `generic` folder, so that you can modify the ones created by the `init.sh` as you please (if you don't remember the syntax, you can go have a look at the ones in the `generic` directory).

##### Installation

To install GODAR simply run the following commands, provided you installed coretran to `GODAR/../coretran`:

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/coretran/lib/cmake
cmake --build build
```

If you want to recompile the code because you made changes, it is recommended that you remove the `build/` directory by `rm -rf build/`.

Something that might streamline the compilation process is to create a file `compile.sh` and add `chmod +x compile.sh` to it. Therefore you can recompile the code with a single command.

```bash
#!/bin/bash
echo 'Removing build files in build/'
rm -rf build/
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/coretran/lib/cmake
cmake --build build
```

#### How to run godar

After compilation, you can use the `start.sh` file if you are working on your own machine. It contains the following line:

```bash
mpirun --bind-to none -n 1 ./godar < input_restart > out
```

The `--bind-to none` makes all threads available to the program (you can then specify this number in the openmp part), while the `-n 1` simply tells how many copies of the program you want to run (on different machines for example). For example, let's say your machine has 2 cores 12 threads per core, you could either run 1 process on 24 threads using openmp or 2 processes (1 per core) with 12 threads using both mpi and openmp.
When working on clusters, you will usually encounter scheduler like slurm. In that case, don't use the `start.sh` file.

The input_file is a simple file to pass along to the main program when executing that can take the following arguments in this order:

- NEED`[bool]   input namelist`
- OPT `[str]    namelist name`
- NEED`[bool]   input restart`
- OPT `[int]    restart version`
- OPT `[int]    restart time`
- NEED`[int]    exp version`
- NEED`[int]    number of particles`
- NEED`[int]    number of threads`

The number of particles must match the number of lines in the initial condition files (or the number of values per line in the output files when restarting). The restart time is the line to restart from in the output files of the restart version.

### How to modify the parameters in GODAR

All physical and numerical parameters as well as options are contained in the namelist.nml file. Therefore, if you want to play with the physics, you can simply change the values in this file without recompiling the code; this permits the use of batch job using an appropriate bash script that modifies the parameters that you want. If you forget the default parameters, they are all set in the `get_default` subroutine in the par_get.f90 file. You can also opt to not use the namelist by setting the read namelist option to false in the input file.

The number of particles is not set at compilation: it is read from the input file, and all the arrays are allocated at run time (`model/src/mod_alloc.f90`), so the same executable works for any number of particles. The only compile-time constants left are in `model/src/mod_parameter.f90`: the size of the legacy rectangular domain (used when no mask is given) and the rank of the master process.

The memory needed per MPI process is about 100n² bytes for n particles, because the state kept for each pair of particles between time steps (contact history and bonds) is stored in n×n arrays. For example, 3 000 particles need about 0.9 GB per process and 10 000 particles about 10 GB, so set the memory of your jobs accordingly.

The input files are setting the run informations: experiment number, whether to use the namelist or not, whether restarting from a previous experiement or not, the number of particles, etc.

### How to run large number of simulations using SLURM

Usually, you will run GODAR on a super computer, which requires a job scheduler like slurm. The scripts in /utils have been written with that usage in mind.

First, you can run `./run_init.py init_args.dat`, which will create initial conditions for your runs. It reads the input_init.dat file, where you can enter various distributions and parameters for both radius and thickness, the rest is handled by the program.

Second, you can run `./duplication.sh < args.dat` which will create namelists based on the arguments provided in the input file. This is for when you want to run the same experiment with multiple types forcings, or different experiments with the same forcing.

Third, you simply run the slurm file the previous program created in `/jobs`. This is where you can modify the allocated time and memory, etc.

Finally, when your jobs have finished running, you can run `./run_video.sh video_args.dat` to create all the videos.

### How to test GODAR

The `tests/` folder contains regression tests, to check that a change to the code does not change the results. After compiling, run them with

```bash
cd build && ctest --output-on-failure
```

or directly with `python tests/run_tests.py` (add `--case TEST` for one case, `-v` to see every comparison). They only need Python 3 (no extra package) and `mpirun`, take about a minute, and need about 2 GB of memory per MPI process.

There are two cases, both with the 6261 floes of the CHANNEL experiment: `TEST` starts from its initial conditions (bonds, winds, sheltering, diagnostics, mask), and `TESTSNAP` from its state at output line 33, when the ice is jammed in the channel (contacts, ridging, walls). Each case runs godar a few times (1 and 2 MPI processes, 1 and 4 OpenMP threads, diagnostics on and off, a restart) and checks that:

- the outputs match the stored references in `tests/reference/`;
- runs that only differ by the number of MPI processes or OpenMP threads give the same results;
- turning the diagnostics off, or restarting, does not change the results at all (bit for bit);
- the ice mass is conserved.

Runs on several processes or threads add up the forces in a different order, which changes the last digits of the results. These comparisons use a relative tolerance of 1e-9 (the differences are about 1e-15 after the 200 steps of a test), small enough to catch any real change of the physics. The references were made on macOS (Apple M3, GNU Fortran 16.2, Open MPI 5.0); other machines and compilers also only change the last digits. The runs and tolerances are set in `tests/cases.json`.

If you change the physics on purpose, the reference check fails: once you are sure the new results are right, store them with `python tests/run_tests.py --update-reference` and commit them with the reason of the change. The run times of each test session are appended to `tests/work/timings.csv`, to follow the performance of the code over time.

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
    - CMakeLists.txt
    - files/
    - generic/
        - input_genreric
        - input_restart_genreric
        - namelist_genreric.nml
        - start_genreric.sh
    - godarCMakesFiles
        - crayFLAGS.cmake
        - fortranENV.cmake
        - gfortranFLAGS.cmake
        - intelFLAGS.cmake
    - model/
        - src/
            - mod_alloc.f90          (allocation of all the module arrays)
            - mod_bond.f90           (module bonds)
            - mod_const.f90          (module const)
            - mod_diagnostics.f90    (module diagnostics)
            - mod_forcings.f90       (module forcings)
            - mod_kind_parameter.f90 (module kind_parameter)
            - mod_mpi.f90            (module mpi_var)
            - mod_options.f90        (module options)
            - mod_pair.f90           (module pairs)
            - mod_parameter.f90      (module parameters)
            - mod_thermo_var.f90     (module thermo_var)
            - mod_variables.f90      (module variables)
            - bonds.f90
            - boundaries.f90
            - compaction.f90
            - contact_forces.f90
            - diagnostics_calc.f90
            - diagnostics_ini_io.f90
            - diagnostics_mpi.f90
            - diagnostics_outputs.f90
            - diagnostics_utils.f90
            - forcing.f90
            - godar.f90
            - ini_get.f90
            - integration.f90
            - mpi_coms.f90
            - par_get.f90
            - physics.f90
            - reset.f90
            - ridging.f90
            - stepper.f90
            - test.f90
            - thermo.f90
    - pkg
        - datetime/
            - date_mod.f90
            - datetime_mod.f90
            - datetimedelta_mod.f90
            - dateutils_mod.f90
            - precision_mod.f90
            - time_mod.f90
        - kdtree/
            - global_kdtree.f90
            - kdtree_utils.f90
        - mask/
            - mask_io.f90
        - mpi/
            - mpi_counts_mod.f90
    - tests/
        - CMakeLists.txt
        - cases.json             (test cases, runs, checks and tolerances)
        - run_tests.py
        - files/                 (initial conditions and mask of the tests)
        - namelist/
        - reference/             (reference outputs of each case)
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
        - init.sh
        - run_init.sh
        - run_video.sh
    - LICENSE
    - README.md
```
