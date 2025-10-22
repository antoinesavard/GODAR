# =================================
# Set cray compile flags
# =================================
message(STATUS "Getting ftn flags")

# Set flags for all build type
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -fbackslash")

if(${LINUX})
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -hfp3 -hvector3 -haggress -hcache3")

  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O0 -g -hflex_mp=conservative -Rb -Rp")
endif()