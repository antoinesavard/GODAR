
message(STATUS "Getting ifort flags")

if(${LINUX})
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -axCORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSSE3 -no-prec-div -fp-model fast=2")

  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -traceback -CB -fp-stack-check -gen-interfaces -warn interfaces")
endif()

if(APPLE)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -axCORE-AVX2,CORE-AVX-I,AVX,SSE4.2,SSSE3 -no-prec-div -fp-model fast=2")

    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -traceback -CB -fp-stack-check -gen-interfaces -warn interfaces")
endif()