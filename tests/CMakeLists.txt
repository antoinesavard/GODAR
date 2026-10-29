# ====================================================
# Regression tests, run after the build with
#     cd build && ctest --output-on-failure
# or directly with python tests/run_tests.py
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
