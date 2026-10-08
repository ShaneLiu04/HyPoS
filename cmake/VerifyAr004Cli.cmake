# VerifyAr004Cli.cmake — AR004 T003 CLI-level verification of the
# --residual-check-interval flag of the hypos binary.
#
# Usage:
#   cmake -Dexe=<path-to-hypos> -Ddir=<work-dir> -Dmode=<mode> \
#         -P VerifyAr004Cli.cmake
#
# Modes:
#   interval10 — jacobi e2e with --residual-check-interval 10, unreachable
#                tolerance, --max-iter 300 on a 32x32 grid: the report JSON
#                must carry "residual_check_interval": 10 and
#                "iterations": N; the profiler output must show
#                residual_allreduce invoked N/10 times (integer division)
#                and residual_confirm invoked exactly once.
#   default    — without the flag the JSON default is 1; the CSV row must
#                have 17 fields with overlapComm at index 10 ("0") and
#                residualCheckInterval at index 11 ("1").
#   cg_warn    — with --solver cg and an interval != 1 rank 0 must WARN
#                that the interval is ignored, and the exit code stays 0.
#   help       — --help must list --residual-check-interval.
#
# Any violation terminates with FATAL_ERROR (non-zero exit).

if(NOT DEFINED exe OR NOT DEFINED dir OR NOT DEFINED mode)
    message(FATAL_ERROR
        "usage: cmake -Dexe=<hypos> -Ddir=<dir> -Dmode=<interval10|default|cg_warn|help> "
        "-P VerifyAr004Cli.cmake")
endif()

if(mode STREQUAL "interval10")
    file(MAKE_DIRECTORY "${dir}")
    execute_process(
        COMMAND "${exe}" --nx 32 --ny 32 --max-iter 300 --tol 1e-8
                --residual-check-interval 10 --enable-profiling
                --output-dir "${dir}"
        WORKING_DIRECTORY "${dir}"
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "ar004 cli verify (interval10): hypos exited ${rc}; output: ${out}; error: ${err}")
    endif()

    if(NOT EXISTS "${dir}/performance_report.json")
        message(FATAL_ERROR "ar004 cli verify (interval10): performance_report.json missing in ${dir}; output: ${out}; error: ${err}")
    endif()
    file(READ "${dir}/performance_report.json" json)
    if(NOT json MATCHES "\"residual_check_interval\": 10")
        message(FATAL_ERROR "ar004 cli verify (interval10): report JSON lacks \"residual_check_interval\": 10; json: ${json}")
    endif()

    string(REGEX MATCH "\"iterations\": ([0-9]+)" _iter_match "${json}")
    if("${CMAKE_MATCH_1}" STREQUAL "")
        message(FATAL_ERROR "ar004 cli verify (interval10): could not extract iterations from report; json: ${json}")
    endif()
    set(iterations "${CMAKE_MATCH_1}")
    math(EXPR expected_allreduce "${iterations} / 10")

    string(REGEX MATCH "\"residual_allreduce\": {[^}]*\"calls\": ([0-9]+)" _ar_match "${out}")
    if("${CMAKE_MATCH_1}" STREQUAL "")
        message(FATAL_ERROR "ar004 cli verify (interval10): residual_allreduce calls not found in profiler output (expected ${expected_allreduce}); output: ${out}")
    endif()
    if(NOT CMAKE_MATCH_1 EQUAL expected_allreduce)
        message(FATAL_ERROR "ar004 cli verify (interval10): residual_allreduce calls=${CMAKE_MATCH_1}, expected ${expected_allreduce} (iterations=${iterations}/10); output: ${out}")
    endif()

    string(REGEX MATCH "\"residual_confirm\": {[^}]*\"calls\": ([0-9]+)" _confirm_match "${out}")
    if("${CMAKE_MATCH_1}" STREQUAL "")
        message(FATAL_ERROR "ar004 cli verify (interval10): residual_confirm calls not found in profiler output (expected 1); output: ${out}")
    endif()
    if(NOT CMAKE_MATCH_1 EQUAL 1)
        message(FATAL_ERROR "ar004 cli verify (interval10): residual_confirm calls=${CMAKE_MATCH_1}, expected 1; output: ${out}")
    endif()

elseif(mode STREQUAL "default")
    execute_process(
        COMMAND "${exe}" --nx 32 --ny 32 --max-iter 5 --output-dir "${dir}"
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "ar004 cli verify (default): hypos (json) exited ${rc}; output: ${out}; error: ${err}")
    endif()
    if(NOT EXISTS "${dir}/performance_report.json")
        message(FATAL_ERROR "ar004 cli verify (default): performance_report.json missing in ${dir}; output: ${out}; error: ${err}")
    endif()
    file(READ "${dir}/performance_report.json" json)
    if(NOT json MATCHES "\"residual_check_interval\": 1")
        message(FATAL_ERROR "ar004 cli verify (default): report JSON lacks default \"residual_check_interval\": 1; json: ${json}")
    endif()

    execute_process(
        COMMAND "${exe}" --nx 32 --ny 32 --max-iter 5 --output-format csv
                --output-dir "${dir}"
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "ar004 cli verify (default): hypos (csv) exited ${rc}; output: ${out}; error: ${err}")
    endif()
    if(NOT EXISTS "${dir}/performance_report.csv")
        message(FATAL_ERROR "ar004 cli verify (default): performance_report.csv missing in ${dir}; output: ${out}; error: ${err}")
    endif()
    file(READ "${dir}/performance_report.csv" csv)
    if("${csv}" STREQUAL "")
        message(FATAL_ERROR "ar004 cli verify (default): performance_report.csv is empty; output: ${out}")
    endif()
    string(REPLACE "\n" ";" csv_lines "${csv}")
    list(GET csv_lines 0 csv_line)
    string(REPLACE "," ";" csv_fields "${csv_line}")
    list(LENGTH csv_fields nfields)
    if(NOT nfields EQUAL 17)
        message(FATAL_ERROR "ar004 cli verify (default): CSV row has ${nfields} fields, expected 17; line: ${csv_line}")
    endif()
    list(GET csv_fields 10 overlap_comm)
    list(GET csv_fields 11 residual_check_interval)
    if(NOT overlap_comm STREQUAL "0")
        message(FATAL_ERROR "ar004 cli verify (default): CSV field 10 (overlap_comm) is \"${overlap_comm}\", expected \"0\"; line: ${csv_line}")
    endif()
    if(NOT residual_check_interval STREQUAL "1")
        message(FATAL_ERROR "ar004 cli verify (default): CSV field 11 (residual_check_interval) is \"${residual_check_interval}\", expected \"1\"; line: ${csv_line}")
    endif()

elseif(mode STREQUAL "cg_warn")
    execute_process(
        COMMAND "${exe}" --nx 32 --ny 32 --solver cg --max-iter 10
                --residual-check-interval 5 --output-dir "${dir}"
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "ar004 cli verify (cg_warn): hypos exited ${rc}; output: ${out}; error: ${err}")
    endif()
    if(NOT out MATCHES "residual-check-interval")
        message(FATAL_ERROR "ar004 cli verify (cg_warn): output lacks \"residual-check-interval\" WARN; output: ${out}")
    endif()
    if(NOT out MATCHES "ignored")
        message(FATAL_ERROR "ar004 cli verify (cg_warn): output lacks \"ignored\" WARN; output: ${out}")
    endif()

elseif(mode STREQUAL "help")
    execute_process(
        COMMAND "${exe}" --help
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT out MATCHES "residual-check-interval")
        message(FATAL_ERROR "ar004 cli verify (help): --help output lacks --residual-check-interval; output: ${out}; error: ${err}")
    endif()

else()
    message(FATAL_ERROR "ar004 cli verify: unknown mode '${mode}' (expected interval10|default|cg_warn|help)")
endif()

message(STATUS "ar004 cli verify (${mode}): ok")
