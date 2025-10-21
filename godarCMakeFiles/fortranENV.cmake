# ================================
# CMAKE Script for setting up the fortran compiling and linking flags for different operating systems and compilers

# Call this cmake file from your project's CMakeLists.txt using 'include(/path/to/fortranENV.cmake)'

# ++++++++++++++++++++++++++++++++


enable_language(Fortran)

# Check if linux
if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/run" CACHE PATH "Installation directory" FORCE)
endif()

# ================================
# Default the Build type to RELEASE
# ================================
# Make sure the build type is uppercase
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

if(BT STREQUAL "RELEASE")
    set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are DEBUG, or RELEASE."
        FORCE)
elseif(BT STREQUAL "DEBUG")
    set(CMAKE_BUILD_TYPE DEBUG CACHE STRING
        "Choose the type of build, options are DEBUG, or RELEASE."
        FORCE)
elseif(NOT BT)
    set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are DEBUG, or RELEASE."
        FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE.")
else()
    message(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, or RELEASE.")
endif(BT STREQUAL "RELEASE")
# ++++++++++++++++++++++++++++++++

# ================================
# Find packages and add compiler-specific flags
# ================================
# openmp
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# mpi
find_package(MPI REQUIRED Fortran)

# coretran
find_package(coretran REQUIRED CONFIG)

# netcdf
# in an ideal world, we would use nf-config --prefix then use
# find_library(NETCDF_NAME NAMES netcdff) but it may be possible that 
# prefix is too high in the tree and that we end up with a netcdff that
# does not match the includedir. On the contrary, it is possible that 
# --flibs is empty, in which case ¯\_(ツ)_/¯
find_program(NF_CONFIG nf-config)
if(NF_CONFIG)
    execute_process(COMMAND ${NF_CONFIG} --includedir
        OUTPUT_VARIABLE NETCDF_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    execute_process(COMMAND ${NF_CONFIG} --flibs
        OUTPUT_VARIABLE NETCDF_FLIBS_STR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        )

    # use --prefix if --flibs is empty
    if (${NETCDF_FLIBS_STR} STREQUAL "")
        message(STATUS "nf-config --flibs is empty! Using --prefix instead")

        execute_process(COMMAND ${NF_CONFIG} --prefix
        OUTPUT_VARIABLE NETCDF_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        set(NETCDF_INCLUDE_DIR "${NETCDF_PREFIX}/include")
        set(NETCDF_LIB_DIR "${NETCDF_PREFIX}/lib")

        find_library(NETCDF_LIBRARY NAMES netcdff
                    HINTS ${NETCDF_LIB_DIR})

    else()
        # Extract -L and -l tokens
        set(NETCDF_LIBRARY "")
        string(REGEX MATCHALL "(-L[^ ]+)" NETCDF_LIBRARY "${NETCDF_FLIBS_STR}")
        string(REGEX MATCHALL "(-l[^ ]+)" NETCDF_LINKS "${NETCDF_FLIBS_STR}")

        # only keep "-lnetcdff"
        set(FILTERED_NETCDF_LINKS "")
        foreach(link IN LISTS NETCDF_LINKS)
            string(STRIP "${link}" link_stripped)
            if(NOT link_stripped STREQUAL "-lnetcdff")
                message(STATUS "Removing link: ${link_stripped}")
            else()
                list(APPEND FILTERED_NETCDF_LINKS "${link_stripped}")
            endif()
        endforeach()

        # combine both lists
        foreach(p IN LISTS FILTERED_NETCDF_LINKS)
            list(APPEND NETCDF_LIBRARY "${p}")
        endforeach()
    endif()
    
    message(STATUS "Found NetCDF-Fortran via nf-config:")
    message(STATUS "  Include:  ${NETCDF_INCLUDE_DIR}")
    message(STATUS "  Linking:  ${NETCDF_LIBRARY}")
else()
    # Fallback to manual search
    find_path(NETCDF_INCLUDE_DIR netcdf.mod
    HINTS ${NETCDF_DIR}/include /usr/local/include /usr/include /opt/homebrew/include /opt /cvmfs
    )
    find_library(NETCDF_LIBRARY NAMES netcdff
        HINTS ${NETCDF_DIR}/lib /usr/local/lib /usr/lib /opt/homebrew/lib /opt /cvmfs
    )
    if(NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
        message(STATUS "Found NetCDF-Fortran:")
        message(STATUS "  Include:  ${NETCDF_INCLUDE_DIR}")
        message(STATUS "  Library:  ${NETCDF_LIBRARY}")
    else()
        message(FATAL_ERROR "Could not find NetCDF Fortran library. Try setting -DNETCDF_DIR=/path/to/netcdf-fortran or use nf-config --prefix to get the path")
        message(STATUS "  Include:  ${NETCDF_INCLUDE_DIR}")
        message(STATUS "  Library:  ${NETCDF_LIBRARY}")
    endif()
endif()

# find_package(HDF5 COMPONENTS Fortran)
# if(HDF5_FOUND)
#     list(APPEND NETCDF_LIBRARY ${HDF5_LIBRARIES})
# endif()

# ++++++++++++++++++++++++++++++++

# ================================
# Set gfortran compile flags
# ================================
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    include(${CMAKE_MODULE_PATH}/godarCMakeFiles/gfortranFLAGS.cmake)

# ================================
# Set INTEL compile flags
# ================================
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    include(${CMAKE_MODULE_PATH}/godarCMakeFiles/intelFLAGS.cmake)

# ================================
# Set Cray compile flags
# ================================
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    include(${CMAKE_MODULE_PATH}/godarCMakeFiles/crayFLAGS.cmake)
endif()

# ================================
# Set the output directories for compiled libraries and module files.
# Paths are relative to the build folder where you call cmake
# ================================
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../lib) # Static library location
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../lib) # Shared library location
# Place module files in specific include folder
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/../include)
# Place executables in bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../bin)
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

# ++++++++++++++++++++++++++++++++

# ================================
# Display information to the user
# ================================
message(STATUS "Build type is ${CMAKE_BUILD_TYPE} use option -DCMAKE_BUILD_TYPE=[DEBUG RELEASE] to switch")

if(BT STREQUAL "RELEASE")
    message(STATUS "Using the following compile flags \n ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif(BT STREQUAL "DEBUG")
    message(STATUS "Using the following compile flags \n ${CMAKE_Fortran_FLAGS_DEBUG}")
else()
    message(STATUS "Using the following compile flags \n ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

message(STATUS "THIS STEP IS NOT REQUIRED")
message(STATUS "cmake --install build would install to ${CMAKE_INSTALL_PREFIX}")
message(STATUS "To change, use option --prefix=/path/to/install/to")