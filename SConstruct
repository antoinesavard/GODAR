from glob import glob

# This file sets up the compilation environment for the sea ice model and its libraries.

# Refer questions to
# David Huard <david.huard@gmail.com>
# McGill U.

# Sept. 18, 2008


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
FC = F90 = "gfortran"

env = Environment(
    LIBPATH=[libs],
    FORTRANMODDIR=[include],
    FORTRANMODDIRPREFIX="-J",
    F90PATH=[include],
    F90FLAGS=[
        "-Wall",
        "-Wno-tabs",
        "-fPIC",
        "-fbackslash",
        "-O3",
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

if int(DEBUG):
    env.Append(FORTRANFLAGS="-g")
    env.Append(FORTRANFLAGS="-static")
    env.Append(FORTRANFLAGS="-fbacktrace")
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
