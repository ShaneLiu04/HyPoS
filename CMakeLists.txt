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
    src/solver/pipelined_cg_solver.cpp
    src/solver/mg_two_level_solver.cpp
    src/solver/mg_hierarchy.cpp
    src/solver/mg_pcg.cpp
    src/comm/p2p_exchanger.cpp
    src/comm/datatype_exchanger.cpp
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
        tests/test_mg_two_level.cpp
        tests/test_mg_vcycle.cpp
        tests/test_pipelined_cg.cpp
        tests/test_performance.cpp
        tests/test_mpiio_output.cpp
    )
    target_link_libraries(test_hypos PRIVATE hypos_core GTest::gtest)

    set(HYPOS_TEST_ENV
        "OMP_NUM_THREADS=4;ASAN_OPTIONS=detect_leaks=0;OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")
    # MPI-launched tests run one OpenMP thread per rank: on small CI runners
    # (4 cores) oversubscribing ranks x threads collapses libgomp regions
    # (measured ~100 ms/iteration vs ~5 us) and times the suite out.
    set(HYPOS_TEST_ENV_MPI
        "OMP_NUM_THREADS=1;ASAN_OPTIONS=detect_leaks=0;OMPI_ALLOW_RUN_AS_ROOT=1;OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1")

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
        # against; MPICH's mpiexec rejects --oversubscribe (OpenMPI's
        # accepts it). CI pins -DMPIEXEC_EXECUTABLE=... per leg, so honor
        # that value first; otherwise fall back to name-based detection.
        set(HYPOS_MPIEXEC "${MPIEXEC_EXECUTABLE}")
        set(HYPOS_MPIEXEC_FLAGS "")
        if(MPI_CXX_COMPILER MATCHES "mpich" OR MPI_CXX_LIBRARIES MATCHES "mpich"
           OR MPIEXEC_EXECUTABLE MATCHES "mpich")
            if(NOT MPIEXEC_EXECUTABLE MATCHES "mpich")
                find_program(HYPOS_MPIEXEC_MPICH NAMES mpiexec.mpich mpirun.mpich)
                if(HYPOS_MPIEXEC_MPICH)
                    set(HYPOS_MPIEXEC "${HYPOS_MPIEXEC_MPICH}")
                else()
                    message(WARNING "MPICH detected but mpiexec.mpich not found; MPI tests may fail")
                endif()
            endif()
        else()
            set(HYPOS_MPIEXEC_FLAGS --oversubscribe)
        endif()
        message(STATUS "HyPoS MPI test launcher: ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS}")

        add_test(NAME halo_2d_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform4RanksHalos:MpiEnvTest.*)
        add_test(NAME halo_3d_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform8RanksHalos3D:MpiEnvTest.*)
        # AR006 (design §4.5): named-filter entries for the datatype and true
        # collective halo paths, mirroring halo_2d_mpi/halo_3d_mpi so the
        # np-hardcoded cases only run at their intended world size.
        add_test(NAME datatype_np4
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform4RanksHalosDatatype:MpiEnvTest.*)
        add_test(NAME datatype_np8
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform8RanksHalos3DDatatype:MpiEnvTest.*)
        add_test(NAME collective_np4
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform4RanksHalosCollective:MpiEnvTest.*)
        add_test(NAME collective_np8
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 8
                    $<TARGET_FILE:test_hypos> --gtest_filter=HaloExchangeTest.MpiNonUniform8RanksHalos3DCollective:MpiEnvTest.*)
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
            ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 60)
        # AR005 (design §4.1 FP1): standalone process-count probe entry.
        # Makes the "multi-rank launch silently degraded to one process"
        # failure mode (mpich 4.2.0 on Ubuntu 24.04) a deterministic FAIL
        # instead of a false green.
        add_test(NAME mpi_size_probe
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=MpiEnvTest.*)
        set_tests_properties(mpi_size_probe PROPERTIES
            ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 60)
        add_test(NAME collective_mpi
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --comm-mode collective --nx 32 --ny 32 --max-iter 20
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_collective)
        # AR007 (design §4.5): mg2 solver paths at their intended world sizes.
        add_test(NAME mg2_np1
            COMMAND $<TARGET_FILE:test_hypos> --gtest_filter=TwoLevelMGUnitTest.*:MpiEnvTest.*)
        add_test(NAME mg2_np4
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=TwoLevelMGMpiTest.*:MpiEnvTest.*)
        add_test(NAME mg2_e2e
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --solver mg2 --nx 32 --ny 32 --max-iter 20
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_mg2)
        # AR008 (design §4.5): mgv/mgcg solver paths at their intended world
        # sizes. mgcg's np1/np4 gtest cases live inside the VCycleMG* filters
        # below (no separate entries — entry-count control, mirroring mg2).
        add_test(NAME mgv_np1
            COMMAND $<TARGET_FILE:test_hypos> --gtest_filter=VCycleMGUnitTest.*:MpiEnvTest.*)
        add_test(NAME mgv_np4
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=VCycleMGMpiTest.*:MpiEnvTest.*)
        add_test(NAME mgv_e2e
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --solver mgv --nx 32 --ny 32 --max-iter 20
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_mgv)
        add_test(NAME mgcg_e2e
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --solver mgcg --nx 32 --ny 32 --max-iter 20
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_mgcg)
        # AR009 (design §4.5/§6): pipelined CG (pcg) + residual history.
        # P-series unit tests live in test_pipelined_cg.cpp behind the
        # PipelinedCG* gtest filters; E-series artifact verifiers are the
        # ResidualHistory* / PipelinedE2E gtest suites reading the fixture
        # directories populated by the residual_history_* entries below.
        add_test(NAME pcg_np1
            COMMAND $<TARGET_FILE:test_hypos> --gtest_filter=PipelinedCGUnitTest.*:MpiEnvTest.*)
        add_test(NAME pcg_np4
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:test_hypos> --gtest_filter=PipelinedCGMpiTest.*:MpiEnvTest.*)
        # E5 (design §6): exit code + "PipelinedCG finished" line + report
        # presence via the script; the float clause (final_residual <=
        # tol*(1+5e-4)) runs in the gtest verifier below.
        add_test(NAME pcg_e2e
            COMMAND ${CMAKE_COMMAND} -Dexe=$<TARGET_FILE:hypos>
                    -Ddir=${CMAKE_BINARY_DIR}/ar009_pcg_e2e -Dmode=pcg_e2e
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr009Cli.cmake)
        add_test(NAME pcg_e2e_verify
            COMMAND $<TARGET_FILE:test_hypos> --gtest_filter=PipelinedE2ETest.*:MpiEnvTest.*)
        # E1/E2 (design §6): per-solver residual-history runs (np1, 32x32,
        # tol 1e-6). jacobi goes through the script so its stdout (the
        # "Converged at iteration" line) is captured for the E2① 6-digit
        # comparison; the other solvers run directly.
        add_test(NAME residual_history_jacobi
            COMMAND ${CMAKE_COMMAND} -Dexe=$<TARGET_FILE:hypos>
                    -Ddir=${CMAKE_BINARY_DIR}/ar009_residual_history/jacobi -Dmode=rh_jacobi
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr009Cli.cmake)
        add_test(NAME residual_history_red_black_gs
            COMMAND $<TARGET_FILE:hypos> --solver red_black_gs --nx 32 --ny 32 --tol 1e-6
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history/red_black_gs
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history/red_black_gs/residual_red_black_gs.csv)
        add_test(NAME residual_history_cg
            COMMAND $<TARGET_FILE:hypos> --solver cg --nx 32 --ny 32 --tol 1e-6
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history/cg
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history/cg/residual_cg.csv)
        add_test(NAME residual_history_mg2
            COMMAND $<TARGET_FILE:hypos> --solver mg2 --nx 32 --ny 32 --tol 1e-6
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history/mg2
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history/mg2/residual_mg2.csv)
        add_test(NAME residual_history_mgv
            COMMAND $<TARGET_FILE:hypos> --solver mgv --nx 32 --ny 32 --tol 1e-6
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history/mgv
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history/mgv/residual_mgv.csv)
        add_test(NAME residual_history_mgcg
            COMMAND $<TARGET_FILE:hypos> --solver mgcg --nx 32 --ny 32 --tol 1e-6
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history/mgcg
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history/mgcg/residual_mgcg.csv)
        add_test(NAME residual_history_pcg
            COMMAND $<TARGET_FILE:hypos> --solver pcg --nx 32 --ny 32 --tol 1e-6
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history/pcg
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history/pcg/residual_pcg.csv)
        # E4 (design §6): interval k=10 row-count semantics; max-iter 100
        # truncation keeps the iteration count deterministic (rows =
        # floor(100/10) = 10).
        add_test(NAME residual_history_interval
            COMMAND $<TARGET_FILE:hypos> --solver jacobi --nx 32 --ny 32 --tol 1e-8 --max-iter 100
                    --residual-check-interval 10
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history/interval
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history/interval/residual_interval.csv)
        # E7 (design §6): --save-interval + --residual-history combo (both
        # branches of the composed single-slot callback); the binary
        # backend makes the k-spaced snapshots land on disk.
        add_test(NAME residual_history_save_combo
            COMMAND $<TARGET_FILE:hypos> --solver jacobi --nx 32 --ny 32 --tol 1e-8 --max-iter 20
                    --output-format binary --save-interval 10
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history/save_combo
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history/save_combo/residual_save_combo.csv)
        add_test(NAME residual_history_e2e_verify
            COMMAND $<TARGET_FILE:test_hypos> --gtest_filter=ResidualHistoryE2ETest.*:MpiEnvTest.*)
        # E3 (design §6): np4 residual history — single rank0 CSV, np1
        # trajectory agreement (the residual_history_cg entry above is the
        # np1 reference: identical loadout).
        add_test(NAME residual_history_np4
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --solver cg --nx 32 --ny 32 --tol 1e-6
                    --output-dir ${CMAKE_BINARY_DIR}/ar009_residual_history_np4
                    --residual-history ${CMAKE_BINARY_DIR}/ar009_residual_history_np4/residual_cg.csv)
        add_test(NAME residual_history_np4_verify
            COMMAND $<TARGET_FILE:test_hypos> --gtest_filter=ResidualHistoryNp4E2ETest.*:MpiEnvTest.*)
        # E6 (design §6): unwritable --residual-history paths (nonexistent
        # parent dir + file-as-parent) — non-zero exit with an [ERROR] line.
        add_test(NAME residual_history_badpath
            COMMAND ${CMAKE_COMMAND} -Dexe=$<TARGET_FILE:hypos>
                    -Ddir=${CMAKE_BINARY_DIR}/ar009_badpath -Dmode=badpath
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr009Cli.cmake)
        # W2 (design §6): --help lists pcg (pipelined wording) and
        # --residual-history.
        add_test(NAME pcg_help_check
            COMMAND ${CMAKE_COMMAND} -Dexe=$<TARGET_FILE:hypos>
                    -Ddir=${CMAKE_BINARY_DIR}/ar009_help -Dmode=help
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr009Cli.cmake)
        # AR006 (design §4.5): e2e through the datatype halo path.
        add_test(NAME datatype_e2e
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
                    $<TARGET_FILE:hypos> --comm-mode datatype --nx 32 --ny 32 --max-iter 20
                    --output-dir ${CMAKE_BINARY_DIR}/ctest_datatype)
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

        # AR005 B1 (design §6.3): np=1 e2e through the vtk and binary
        # backends. The verify steps check the exact file set plus the
        # on-disk payload contracts (appended length header / shard header
        # fields), not just file existence.
        add_test(NAME vtk_e2e
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
                    $<TARGET_FILE:hypos> --output-format vtk --nx 16 --ny 16 --max-iter 20
                    --save-interval 10 --output-dir ${CMAKE_BINARY_DIR}/ctest_vtk)
        add_test(NAME vtk_e2e_verify
            COMMAND ${CMAKE_COMMAND}
                    -Ddir=${CMAKE_BINARY_DIR}/ctest_vtk
                    "-Dexpected=solution_10.pvti;solution_10_r0.vti;solution_20.pvti;solution_20_r0.vti"
                    -Dfmt=vtk -Dpiece=solution_10_r0.vti -Dnx=16 -Dny=16 -Dnz=1
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr005Output.cmake)
        add_test(NAME binary_e2e
            COMMAND ${HYPOS_MPIEXEC} ${HYPOS_MPIEXEC_FLAGS} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 1
                    $<TARGET_FILE:hypos> --output-format binary --nx 16 --ny 16 --max-iter 20
                    --save-interval 10 --output-dir ${CMAKE_BINARY_DIR}/ctest_binary)
        add_test(NAME binary_e2e_verify
            COMMAND ${CMAKE_COMMAND}
                    -Ddir=${CMAKE_BINARY_DIR}/ctest_binary
                    "-Dexpected=solution_10_r0.bin;solution_20_r0.bin"
                    -Dfmt=binary -Dpiece=solution_10_r0.bin -Dnx=16 -Dny=16 -Dnz=1
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifyAr005Output.cmake)

        # AR005 (design §4.3): inject the expected rank count into every
        # np>1 entry. gtest entries consume it via the MpiEnvTest probe
        # (appended to their --gtest_filter above); the four hypos-binary
        # entries are declarative self-documentation only -- their guard is
        # the environment-level probe (decision D7: launch degradation is an
        # all-or-nothing property of the mpirun environment).
        # ENV base is HYPOS_TEST_ENV_MPI (CI fix, commit 4ce4bea): MPI-launched
        # tests run one OpenMP thread per rank to avoid oversubscribing
        # small CI runners.
        set_tests_properties(halo_2d_mpi
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 600)
        set_tests_properties(halo_3d_mpi
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=8" TIMEOUT 600)
        set_tests_properties(datatype_np4 collective_np4
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 600)
        set_tests_properties(datatype_np8 collective_np8
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=8" TIMEOUT 600)
        # Heavy multi-iteration suites get extra headroom for the slow
        # Debug+ASan CI leg; they finish in seconds on normal hardware.
        set_tests_properties(solver_mpi overlap_consistency solver_alt_mpi
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 900)
        set_tests_properties(collective_mpi datatype_e2e
            mpiio_e2e mpiio_e2e_final
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 600)
        set_tests_properties(mg2_np1
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 900)
        set_tests_properties(mg2_np4 mg2_e2e
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 900)
        # AR008 (design §4.5): same ENV/TIMEOUT pattern as the mg2 entries.
        set_tests_properties(mgv_np1
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 900)
        set_tests_properties(mgv_np4 mgv_e2e mgcg_e2e
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 900)
        # AR009 (design §4.5): same ENV/TIMEOUT pattern as the mg2/mgv
        # entries. The e2e run entries FIXTURES_SETUP the artifact
        # directories consumed by the *_verify gtest entries (mpiio
        # fixture pattern); the verify entries carry the artifact path via
        # HYPOS_AR009_* env vars and SKIP when those are unset (plain
        # `unit` entry).
        set_tests_properties(pcg_np1
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 900)
        set_tests_properties(pcg_np4
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 900)
        set_tests_properties(pcg_e2e residual_history_jacobi
            residual_history_badpath pcg_help_check
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 300
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
        set_tests_properties(pcg_e2e
            PROPERTIES FIXTURES_SETUP ar009_pcg_e2e_files)
        set_tests_properties(pcg_e2e_verify
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV};HYPOS_AR009_PCG_E2E_DIR=${CMAKE_BINARY_DIR}/ar009_pcg_e2e"
                FIXTURES_REQUIRED ar009_pcg_e2e_files TIMEOUT 60)
        set_tests_properties(residual_history_jacobi
            PROPERTIES FIXTURES_SETUP ar009_rh_jacobi_files)
        set_tests_properties(residual_history_red_black_gs
            PROPERTIES FIXTURES_SETUP ar009_rh_red_black_gs_files)
        set_tests_properties(residual_history_cg
            PROPERTIES FIXTURES_SETUP ar009_rh_cg_files)
        set_tests_properties(residual_history_mg2
            PROPERTIES FIXTURES_SETUP ar009_rh_mg2_files)
        set_tests_properties(residual_history_mgv
            PROPERTIES FIXTURES_SETUP ar009_rh_mgv_files)
        set_tests_properties(residual_history_mgcg
            PROPERTIES FIXTURES_SETUP ar009_rh_mgcg_files)
        set_tests_properties(residual_history_pcg
            PROPERTIES FIXTURES_SETUP ar009_rh_pcg_files)
        set_tests_properties(residual_history_interval
            PROPERTIES FIXTURES_SETUP ar009_rh_interval_files)
        set_tests_properties(residual_history_save_combo
            PROPERTIES FIXTURES_SETUP ar009_rh_save_combo_files)
        set_tests_properties(residual_history_red_black_gs residual_history_cg
            residual_history_mg2 residual_history_mgv residual_history_mgcg
            residual_history_pcg residual_history_interval residual_history_save_combo
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 300)
        set_tests_properties(residual_history_e2e_verify
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV};HYPOS_AR009_RH_BASE=${CMAKE_BINARY_DIR}/ar009_residual_history"
                FIXTURES_REQUIRED "ar009_rh_jacobi_files;ar009_rh_red_black_gs_files;ar009_rh_cg_files;ar009_rh_mg2_files;ar009_rh_mgv_files;ar009_rh_mgcg_files;ar009_rh_pcg_files;ar009_rh_interval_files;ar009_rh_save_combo_files"
                TIMEOUT 120)
        set_tests_properties(residual_history_np4
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 300
                FIXTURES_SETUP ar009_rh_np4_files)
        set_tests_properties(residual_history_np4_verify
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV};HYPOS_AR009_RH_NP4_DIR=${CMAKE_BINARY_DIR}/ar009_residual_history_np4;HYPOS_AR009_RH_BASE=${CMAKE_BINARY_DIR}/ar009_residual_history"
                FIXTURES_REQUIRED "ar009_rh_np4_files;ar009_rh_cg_files" TIMEOUT 60)
        set_tests_properties(mpiio_np4
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=4" TIMEOUT 600)
        set_tests_properties(ranks3_smoke
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI};HYPOS_EXPECT_NP=3" TIMEOUT 600)
        set_tests_properties(mpiio_np1
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV_MPI}" TIMEOUT 600)
        set_tests_properties(mpiio_e2e PROPERTIES FIXTURES_SETUP mpiio_e2e_files)
        set_tests_properties(mpiio_e2e_verify PROPERTIES FIXTURES_REQUIRED mpiio_e2e_files)
        set_tests_properties(mpiio_e2e_final PROPERTIES FIXTURES_SETUP mpiio_e2e_final_files)
        set_tests_properties(mpiio_e2e_final_verify PROPERTIES FIXTURES_REQUIRED mpiio_e2e_final_files)
        set_tests_properties(vtk_e2e binary_e2e
            PROPERTIES ENVIRONMENT "${HYPOS_TEST_ENV}" TIMEOUT 600)
        set_tests_properties(vtk_e2e PROPERTIES FIXTURES_SETUP vtk_e2e_files)
        set_tests_properties(vtk_e2e_verify PROPERTIES FIXTURES_REQUIRED vtk_e2e_files)
        set_tests_properties(binary_e2e PROPERTIES FIXTURES_SETUP binary_e2e_files)
        set_tests_properties(binary_e2e_verify PROPERTIES FIXTURES_REQUIRED binary_e2e_files)
    endif()
elseif(HYPOS_BUILD_TESTS)
    message(WARNING "GTest not found: tests are disabled")
endif()

# Installation
install(TARGETS hypos DESTINATION bin)
