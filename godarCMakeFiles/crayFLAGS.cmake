
message(STATUS "Getting ftn flags")

if(${LINUX})
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fPIC -fbackslash -O3 -hfp3 -hvector3 -haggress -hcache3 -hconcurrent")

  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Og -g -static -fbacktrace -fcheck=all")
endif()