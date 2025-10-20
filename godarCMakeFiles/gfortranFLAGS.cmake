# ================================
# Set gfortran compile flags
# ================================
message(STATUS "Getting gfortran flags")

# Set flags for all build types
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wno-tabs -fPIC -fbackslash")

set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -funroll-all-loops -finline-functions -ffast-math")

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -Og -fbacktrace -static -fcheck=all")