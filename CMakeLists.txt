cmake_minimum_required(VERSION 3.16)
project(HyPoS VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build types
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Compiler options (GNU / Clang / Intel; explicit failure otherwise)
include(CompilerWarnings)

# OpenMP
find_package(OpenMP)

# MPI
find_package(MPI REQUIRED)

# Source files
set(HYPOS_CORE_SOURCES
    src/grid/subgrid.cpp
    src/grid/partition.cpp
    src/solver/jacobi_solver.cpp
    src/solver/red_black_gs_solver.cpp
    src/solver/cg_solver.cpp
    src/comm/p2p_exchanger.cpp
    src/comm/collective_exchanger.cpp
    src/perf/timer.cpp
    src/perf/profiler.cpp
    src/perf/reporter.cpp
    src/utils/mpi_env.cpp
    src/utils/logger.cpp
    src/utils/cmdline_parser.cpp
    src/io/binary_io.cpp
    src/io/vtk_io.cpp
    src/io/mpiio_binary.cpp
)

# Core library shared by the executable and tests
add_library(hypos_core STATIC ${HYPOS_CORE_SOURCES})
target_include_directories(hypos_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(hypos_core PUBLIC MPI::MPI_CXX)
if(OpenMP_CXX_FOUND)
    target_link_libraries(hypos_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# Main executable
add_executable(hypos src/main.cpp)
target_link_libraries(hypos PRIVATE hypos_core)

# Tests
option(HYPOS_BUILD_TESTS "Build tests when GTest is available" ON)
enable_testing()
if(HYPOS_BUILD_TESTS)
    find_package(GTest)
else()
    set(GTest_FOUND FALSE)
endif()
if(GTest_FOUND)
    add_executable(test_hypos
        tests/gtest_main_mpi.cpp
        tests/test_grid.cpp
        tests/test_solver.cpp
        tests/test_alt_solvers.cpp
        tests/test_overlap.cpp
        tests/test_halo_exchange.cpp
        tests/test_solver_mpi.cpp
        tests/test_performance.cpp
        tests/test_mpiio_output.cpp
    )
    target_link_libraries(test_hypos PRIVATE hypos_core GTest::gtest)

    set(HYPOS_TEST_ENV
        "OMP_NUM_THREADS=4;ASAN_OPTIONS=detect_leaks=0;OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")

    add_test(NAME unit COMMAND $<TARGET_FILE:test_hypos>)
    set_tests_properties(unit PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 120)

    add_test(NAME parse_error COMMAND $<TARGET_FILE:hypos> --save-interval abc)
    set_tests_properties(parse_error PROPERTIES WILL_FAIL TRUE TIMEOUT 60)

    if(MPIEXEC_EXECUTABLE)
        # Launcher must match the MPI implementation the binary is linked
        # against; MPICH's mpiexec does not accept --oversubscribe (which
        # keeps the np=8 test runnable on machines with fewer cores).
        set(HYPOS_MPIEXEC "${MPIEXEC_EXECUTABLE}")
        set(HYPOS_MPIEXEC_FLAGS "")
        if(MPI_CXX_COMPILER MATCHES "mpich")
            find_program(HYPOS_MPIEXEC_MPICH NAMES mpiexec.mpich mpirun.mpich)
            if(HYPOS_MPIEXEC_MPICH)
                set(HYPOS_MPIEXEC "${HYPOS_MPIEXEC_MPICH}")
            endif()
        else()
            set(HYPOS_MPIEXEC_FLAGS --oversubscribe)
        endif()

        add_test(NAME halo_2d_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform4RanksHalos)
        add_test(NAME halo_3d_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform8RanksHalos3D)
        add_test(NAME solver_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=SolverMpiTest.*)
        add_test(NAME overlap_consistency
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=OverlapTest.*)
        add_test(NAME solver_alt_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=SolverMpiTest.AltSolvers*)
        add_test(NAME collective_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --comm-mode collective --nx 32 --ny 32 --max-iter 20
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_collective)
        add_test(NAME ranks3_smoke
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 3
                    $<TARGET_FILE:hypos> --nx 24 --ny 24 --max-iter 20
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_ranks3)
        add_test(NAME mpiio_np1
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
                    $<TARGET_FILE:test_hypos> --gtest_filter=MpiIoOutputTest.*)
        add_test(NAME mpiio_np4
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=MpiIoOutputTest.*)
        add_test(NAME mpiio_e2e
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --output-format mpibin --nx 32 --ny 32 --max-iter 20
                    --save-interval 10 --enable-profiling
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_mpiio)

        set_tests_properties(halo_2d_mpi halo_3d_mpi solver_mpi overlap_consistency
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 600)
        set_tests_properties(solver_alt_mpi collective_mpi ranks3_smoke
            mpiio_np1 mpiio_np4 mpiio_e2e
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 600)
    endif()
elseif(HYPOS_BUILD_TESTS)
    message(WARNING "GTest not found: tests are disabled")
endif()

# Installation
install(TARGETS hypos DESTINATION bin)
