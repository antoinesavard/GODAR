from glob import glob
import os

# This file sets up the compilation environment for Godar and its libraries.

# ----------------------------------
#         SCONS Variables
# ----------------------------------

# LIBPATH : The path searched for libraries.
# FORTRANMODDIR : The path where .mod files are stored, ie compiled modules.
# FORTRANMODDIRPREFIX : The compiler prefix that prefixes FORTRANMODDIR (-J for gfortran)
# F90PATH : The path searched for include files.
# F90FLAGS : Flags for the F90 compiler.


# ----------------------------------
#              Notes
# ----------------------------------
#
# Flags
# -----
#   -fbackslash : Allow escape characters \n, \t, etc in strings.
#                 In gcc-4.3, -fno-backslash is the default.
#   -Wall : Print all warnings.

# ----------------------------------
#  Standard Compilation Environment
# ----------------------------------

include = Dir("#/include")
libs = Dir("#/libs")
FC = F90 = "mpifort"

env = Environment(
    ENV=os.environ,
    LIBPATH=[libs],
    FORTRANMODDIR=[include],
    FORTRANMODDIRPREFIX="-J",
    F90PATH=[include],
    F90FLAGS=[
        "-Wall",
        "-Wno-tabs",
        "-fPIC",
        "-fbackslash",
        "-ffast-math",
        "-mcmodel=medium",
        "-fopenmp",
    ],
    LINKFLAGS=[
        "-fopenmp",
    ],
    FORTRAN=FC,
    F90=FC,
)

env.Append(LIBPATH="/aos/home/asavard/coretran/lib")
env.Append(F90PATH="/aos/home/asavard/coretran/include")

# Command line options to modify the build environment.
# Use scons-3 debug=1 to compile with debugging flags.


DEBUG = ARGUMENTS.get("debug", 0)
EXE = ARGUMENTS.get("exe", "godar")

if int(DEBUG) == 0:
    env.Append(F90FLAGS="-O3")

if int(DEBUG):
    env.Append(F90FLAGS="-g")
    env.Append(F90FLAGS="-Og")
    env.Append(F90FLAGS="-static")
    env.Append(F90FLAGS="-fbacktrace")
    env.Append(F90FLAGS="-fcheck=all")
    env.Append(LINKFLAGS="-fbacktrace")

# ----------------------------------------
#             Create Executables
# ----------------------------------------

# Export the compilation environment so SConscript files in subdirectories can access it.

Export("env")
Export("EXE")

# The datetime module
datetime = SConscript("datetime/SConscript", variant_dir="build/datetime")

# The GoDAR library and executable
godar = SConscript("src/SConscript", variant_dir="build/godar")
