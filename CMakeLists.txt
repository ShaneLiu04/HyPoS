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

# Compiler options
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-march=native>
    $<$<CONFIG:Release>:-DNDEBUG>
    $<$<CONFIG:Debug>:-O0>
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Debug>:-fsanitize=address,undefined>
)

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions(-D_OPENMP)
endif()

# MPI
find_package(MPI REQUIRED)

# PAPI (optional)
find_package(PAPI)
if(PAPI_FOUND)
    add_definitions(-DHYPOS_ENABLE_PAPI)
endif()

# Source files
set(HYPOS_SOURCES
    src/main.cpp
    src/core/memory_pool.cpp
    src/grid/subgrid.cpp
    src/grid/partition.cpp
    src/solver/jacobi_solver.cpp
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
)

# Main executable
add_executable(hypos ${HYPOS_SOURCES})
target_include_directories(hypos PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(hypos PRIVATE MPI::MPI_CXX)
if(OpenMP_CXX_FOUND)
    target_link_libraries(hypos PRIVATE OpenMP::OpenMP_CXX)
endif()
if(PAPI_FOUND)
    target_include_directories(hypos PRIVATE ${PAPI_INCLUDE_DIRS})
    target_link_libraries(hypos PRIVATE ${PAPI_LIBRARIES})
endif()

# Tests
enable_testing()
find_package(GTest)
if(GTest_FOUND)
    add_executable(test_hypos
        tests/test_grid.cpp
        tests/test_solver.cpp
        tests/test_halo_exchange.cpp
        tests/test_performance.cpp
    )
    target_include_directories(test_hypos PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(test_hypos PRIVATE GTest::gtest GTest::gtest_main MPI::MPI_CXX)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(test_hypos PRIVATE OpenMP::OpenMP_CXX)
    endif()
    add_test(NAME HyPoS_UnitTests COMMAND test_hypos)
endif()

# Installation
install(TARGETS hypos DESTINATION bin)
