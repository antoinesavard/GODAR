# ================================
# Set gfortran compile flags
# ================================
message(STATUS "Getting gfortran flags")

# Set flags for all build types
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wno-tabs -fPIC -fbackslash")

set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -funroll-all-loops -finline-functions -ffast-math")

# on APPLE, it is possible that -static doesn't work
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -g -Og -static -fbacktrace -fcheck=all")