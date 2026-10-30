# =================================
# Set AMD Flang (LLVMFlang) compile flags
# =================================
message(STATUS "Getting amdflang flags")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -funroll-loops -ffast-math")

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O0 -g")
