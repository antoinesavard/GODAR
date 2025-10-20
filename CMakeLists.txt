cmake_minimum_required(VERSION 3.8)

project(godar LANGUAGES Fortran)

# Optional parameters to use
# -DCMAKE_PREFIX_PATH=/path/to/coretran
# -DCMAKE_Fortran_COMPILER=/path/to/fortran/compiler
# -DCMAKE_BUILD_TYPE=[DEBUG, RELEASE]  - Build with debug flags or not.
# -DBUILD_SHARED_LIBS=[ON, OFF]  - Build and shared or static library, or link a program using the appropriate flags

# ====================================================
# Setup CMake paths and Fortran environment
# ====================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
include("${CMAKE_CURRENT_SOURCE_DIR}/godarCMakeFiles/fortranENV.cmake")
include(GNUInstallDirs)

# ====================================================
# Build all pkg submodules automatically
# ====================================================
set(PKG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pkg)
file(GLOB PKG_SUBDIRS LIST_DIRECTORIES true ${PKG_DIR}/*)

set(PKG_LIBS "")  # will hold all module libraries

foreach(subdir ${PKG_SUBDIRS})
    if (IS_DIRECTORY ${subdir})
        get_filename_component(modname ${subdir} NAME)
        message(STATUS "Adding module library: ${modname}")

        file(GLOB MOD_SRC ${subdir}/*.f90)
        if (MOD_SRC)
            add_library(${modname} ${MOD_SRC})
            # Include paths for Fortran modules
            target_include_directories(${modname} PUBLIC
                ${CMAKE_CURRENT_SOURCE_DIR}/model/inc
                ${CMAKE_Fortran_MODULE_DIRECTORY}
                ${subdir}
                ${MPI_Fortran_INCLUDE_DIRS}
                ${NETCDF_INCLUDE_DIR}
                ${CORETRAN_INCLUDE_DIR}
            )
            # Link each library with Coretran (and optionally MPI/NetCDF if needed)
            target_link_libraries(${modname} PUBLIC 
                coretran
                ${OpenMP_Fortran_LIBRARIES}
                ${NETCDF_LIBRARY}
                ${MPI_Fortran_LIBRARIES}
            )
            list(APPEND PKG_LIBS ${modname})
        endif()
    endif()
endforeach()

# ====================================================
# Gather all model source files
# ====================================================
file(GLOB MODEL_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/model/src/*.f90)
add_executable(${PROJECT_NAME} ${MODEL_SRC_FILES})

# Include directories for Fortran INCLUDE statements (.h)
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/model/inc
    ${MPI_Fortran_INCLUDE_DIRS}
    ${NETCDF_INCLUDE_DIR}
    ${CORETRAN_INCLUDE_DIR}
)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE  
    coretran
    ${PKG_LIBS}
    ${OpenMP_Fortran_LIBRARIES}
    ${NETCDF_LIBRARY}
    ${MPI_Fortran_LIBRARIES}
)

# ====================================================
# Install targets
# ====================================================
install(TARGETS godar
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/godar
    FILES_MATCHING PATTERN "*.mod"
)