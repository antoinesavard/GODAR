# ====================================================
# Regression tests, run after the build with
#     cd build && ctest --output-on-failure
# ====================================================
find_program(GODAR_PYTHON NAMES python3 python)
if(NOT GODAR_PYTHON)
    message(WARNING "python not found: the regression tests are not added")
    return()
endif()

# test cases of cases.json
set(GODAR_TEST_CASES TEST TESTSNAP)

# MPI launcher found by find_package(MPI) in fortranENV.cmake
set(GODAR_MPIEXEC_ARG "")
if(MPIEXEC_EXECUTABLE)
    set(GODAR_MPIEXEC_ARG --mpiexec ${MPIEXEC_EXECUTABLE})
endif()

foreach(case IN LISTS GODAR_TEST_CASES)
    add_test(NAME ${case}
        COMMAND ${GODAR_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.py
                --godar $<TARGET_FILE:${PROJECT_NAME}>
                ${GODAR_MPIEXEC_ARG}
                --case ${case})
    # each case uses several GB of memory and all the cores: one at a time
    set_tests_properties(${case} PROPERTIES RESOURCE_LOCK godar TIMEOUT 3600)
endforeach()

# benchmark of all the cases (about 5 minutes), compared with the stored
# baseline of the machine if there is one; only run with
#     cd build && ctest -C bench -R bench --verbose
add_test(NAME bench CONFIGURATIONS bench
    COMMAND ${GODAR_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/run_bench.py
            --godar $<TARGET_FILE:${PROJECT_NAME}>
            ${GODAR_MPIEXEC_ARG})
set_tests_properties(bench PROPERTIES RESOURCE_LOCK godar TIMEOUT 7200)

# ARCTIC benchmark: GODAR_BENCH_PARTICLES floes filling the Arctic sea ice
# extent, on 1, 2, 4 nodes x 1, 2, 4 MPI processes per node x 16 to 128
# OpenMP threads; the configurations that do not fit the Slurm allocation
# are skipped. arctic_case builds the case first (kept, only rebuilt if
# missing). Run with
#     cd build && ctest -C bench -R arctic -V
set(GODAR_BENCH_PARTICLES 20000 CACHE STRING
    "number of floes of the ARCTIC benchmark")
set(GODAR_ARCTIC_CASE
    ${CMAKE_CURRENT_SOURCE_DIR}/work/bench/cases/ARCTIC_${GODAR_BENCH_PARTICLES})
add_test(NAME arctic_case CONFIGURATIONS bench
    COMMAND ${GODAR_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/bench/make_arctic.py
            --n ${GODAR_BENCH_PARTICLES} --out ${GODAR_ARCTIC_CASE})
add_test(NAME arctic CONFIGURATIONS bench
    COMMAND ${GODAR_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/run_bench.py
            --godar $<TARGET_FILE:${PROJECT_NAME}>
            ${GODAR_MPIEXEC_ARG}
            --particles ${GODAR_BENCH_PARTICLES}
            --nodes 1,2,4 --ranks-per-node 1,2,4 --threads 16,32,64,128)
set_tests_properties(arctic_case PROPERTIES FIXTURES_SETUP arctic TIMEOUT 3600)
set_tests_properties(arctic PROPERTIES FIXTURES_REQUIRED arctic
                     RESOURCE_LOCK godar TIMEOUT 14400)
