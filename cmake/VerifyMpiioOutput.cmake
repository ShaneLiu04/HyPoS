# VerifyMpiioOutput.cmake — assert that an mpibin output directory contains
# EXACTLY the expected file set (single-file products only: no per-rank
# "_r<rank>" shards, no parallel index, no leftovers from earlier runs).
#
# Usage:
#   cmake -Ddir=<output-dir> -Dexpected="solution_10.bin;solution_20.bin" \
#         -P VerifyMpiioOutput.cmake
#
# Fails (non-zero exit) when the directory is missing, or when the actual
# file list differs from the expected list in any way.

if(NOT DEFINED dir OR NOT DEFINED expected)
    message(FATAL_ERROR "usage: cmake -Ddir=<dir> -Dexpected=<f1;f2> -P VerifyMpiioOutput.cmake")
endif()

if(NOT EXISTS "${dir}")
    message(FATAL_ERROR "output dir does not exist: ${dir}")
endif()

file(GLOB actual_files "${dir}/*")
set(actual_names "")
foreach(f IN LISTS actual_files)
    get_filename_component(fname "${f}" NAME)
    list(APPEND actual_names "${fname}")
endforeach()
list(SORT actual_names)

set(expected_names ${expected})
list(SORT expected_names)

if(NOT "${actual_names}" STREQUAL "${expected_names}")
    message(FATAL_ERROR
        "unexpected mpibin output files in ${dir}\n"
        "  actual:   ${actual_names}\n"
        "  expected: ${expected_names}\n"
        "(single-file format must produce exactly the expected files; "
        "per-rank shards or index files must not appear)")
endif()

message(STATUS "mpibin output verified (${dir}): ${actual_names}")
