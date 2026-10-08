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
    src/solver/residual.cpp
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

    # E-3 (design §6.2/§6.7): an unknown --output-format value must be
    # rejected with a non-zero exit code.
    add_test(NAME mpibin_bad_format COMMAND $<TARGET_FILE:hypos> --output-format mpibinn)
    set_tests_properties(mpibin_bad_format PROPERTIES WILL_FAIL TRUE TIMEOUT 60)

    # AR004 T003: invalid --residual-check-interval usage must be rejected
    # with exit code 1 (WILL_FAIL, mirroring parse_error/mpibin_bad_format).
    add_test(NAME ar004_interval_zero COMMAND $<TARGET_FILE:hypos> --residual-check-interval=0)
    add_test(NAME ar004_interval_neg COMMAND $<TARGET_FILE:hypos> --residual-check-interval=-5)
    add_test(NAME ar004_interval_nan COMMAND $<TARGET_FILE:hypos> --residual-check-interval=abc)
    add_test(NAME ar004_interval_missing COMMAND $<TARGET_FILE:hypos> --residual-check-interval)
    set_tests_properties(ar004_interval_zero ar004_interval_neg ar004_interval_nan ar004_interval_missing
        PROPERTIES WILL_FAIL TRUE TIMEOUT 60)

    # AR004 T003: CLI e2e verification (report JSON field, CSV column
    # position, cg WARN, --help listing) via the script-mode verifier.
    add_test(NAME ar004_cli_interval10
        COMMAND ${CMAKE_COMMAND} -Dexe=$<TARGET_FILE:hypos>
                -Ddir=${CMAKE_BINARY_DIR}/ar004_cli_interval10 -Dmode=interval10
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr004Cli.cmake)
    add_test(NAME ar004_cli_default
        COMMAND ${CMAKE_COMMAND} -Dexe=$<TARGET_FILE:hypos>
                -Ddir=${CMAKE_BINARY_DIR}/ar004_cli_default -Dmode=default
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr004Cli.cmake)
    add_test(NAME ar004_cli_cg_warn
        COMMAND ${CMAKE_COMMAND} -Dexe=$<TARGET_FILE:hypos>
                -Ddir=${CMAKE_BINARY_DIR}/ar004_cli_cg_warn -Dmode=cg_warn
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr004Cli.cmake)
    add_test(NAME ar004_cli_help
        COMMAND ${CMAKE_COMMAND} -Dexe=$<TARGET_FILE:hypos>
                -Ddir=${CMAKE_BINARY_DIR}/ar004_cli_help -Dmode=help
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr004Cli.cmake)
    set_tests_properties(ar004_cli_interval10 ar004_cli_default ar004_cli_cg_warn ar004_cli_help
        PROPERTIES TIMEOUT 60 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        ENVIRONMENT "${HYPOS_TEST_ENV}")

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
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform4RanksHalos:MpiEnvTest.*)
        add_test(NAME halo_3d_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform8RanksHalos3D:MpiEnvTest.*)
        add_test(NAME solver_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=SolverMpiTest.*:MpiEnvTest.*)
        add_test(NAME overlap_consistency
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=OverlapTest.*:MpiEnvTest.*)
        add_test(NAME solver_alt_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=SolverMpiTest.AltSolvers*:MpiEnvTest.*)
        # B3 (AR004 T006, design §6): np=4 topology consistency of the
        # Cartesian communicator exposed by UniformPartition::partition().
        add_test(NAME topology_consistency
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=PartitionTopologyTest.*:MpiEnvTest.*)
        set_tests_properties(topology_consistency PROPERTIES
            ENVIRONMENT "${HYPOS_TEST_ENV};HYPOS_EXPECT_NP=4" TIMEOUT 60)
        # AR005 (design §4.1 FP1): standalone process-count probe entry.
        # Makes the "multi-rank launch silently degraded to one process"
        # failure mode (mpich 4.2.0 on Ubuntu 24.04) a deterministic FAIL
        # instead of a false green.
        add_test(NAME mpi_size_probe
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=MpiEnvTest.*)
        set_tests_properties(mpi_size_probe PROPERTIES
            ENVIRONMENT "${HYPOS_TEST_ENV};HYPOS_EXPECT_NP=4" TIMEOUT 60)
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
                    $<TARGET_FILE:test_hypos> --gtest_filter=MpiIoOutputTest.*:MpiEnvTest.*)
        add_test(NAME mpiio_e2e
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --output-format mpibin --nx 32 --ny 32 --max-iter 20
                    --save-interval 10 --enable-profiling
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_mpiio)
        # S-4 (design §6.3): assert the save-interval run produced exactly
        # the expected single-file snapshots (no shards, no index).
        add_test(NAME mpiio_e2e_verify
            COMMAND ${CMAKE_COMMAND}
                    -Ddir=${CMAKE_BINARY_DIR}/ctest_mpiio
                    "-Dexpected=solution_10.bin;solution_20.bin"
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyMpiioOutput.cmake)
        # S-5 (design §6.3): without --save-interval only the final solution
        # file is written.
        add_test(NAME mpiio_e2e_final
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --output-format mpibin --nx 16 --ny 16 --max-iter 5
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_mpiio_final)
        add_test(NAME mpiio_e2e_final_verify
            COMMAND ${CMAKE_COMMAND}
                    -Ddir=${CMAKE_BINARY_DIR}/ctest_mpiio_final
                    "-Dexpected=solution_5.bin"
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyMpiioOutput.cmake)

        # AR005 (design §4.3): inject the expected rank count into every
        # np>1 entry. gtest entries consume it via the MpiEnvTest probe
        # (appended to their --gtest_filter above); the four hypos-binary
        # entries are declarative self-documentation only -- their guard is
        # the environment-level probe (decision D7: launch degradation is an
        # all-or-nothing property of the mpirun environment).
        set_tests_properties(halo_2d_mpi solver_mpi overlap_consistency
            solver_alt_mpi mpiio_np4 collective_mpi mpiio_e2e mpiio_e2e_final
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV};HYPOS_EXPECT_NP=4" TIMEOUT 600)
        set_tests_properties(halo_3d_mpi
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV};HYPOS_EXPECT_NP=8" TIMEOUT 600)
        set_tests_properties(ranks3_smoke
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV};HYPOS_EXPECT_NP=3" TIMEOUT 600)
        set_tests_properties(mpiio_np1
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 600)
        set_tests_properties(mpiio_e2e PROPERTIES FIXTURES_SETUP mpiio_e2e_files)
        set_tests_properties(mpiio_e2e_verify PROPERTIES FIXTURES_REQUIRED mpiio_e2e_files)
        set_tests_properties(mpiio_e2e_final PROPERTIES FIXTURES_SETUP mpiio_e2e_final_files)
        set_tests_properties(mpiio_e2e_final_verify PROPERTIES FIXTURES_REQUIRED mpiio_e2e_final_files)
    endif()
elseif(HYPOS_BUILD_TESTS)
    message(WARNING "GTest not found: tests are disabled")
endif()

# Installation
install(TARGETS hypos DESTINATION bin)
