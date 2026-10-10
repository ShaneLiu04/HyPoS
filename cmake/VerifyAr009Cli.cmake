# VerifyAr009Cli.cmake — AR009 T001 CLI-level verification of the hypos
# binary (pipelined CG + residual history).
#
# Usage:
#   cmake -Dexe=<path-to-hypos> -Ddir=<work-dir> -Dmode=<mode> \
#         -P VerifyAr009Cli.cmake
#
# Modes:
#   pcg_e2e   — run `--solver pcg` (np1, 32x32, tol 1e-6): exit code 0, the
#               "PipelinedCG finished" log line on stdout, and a
#               performance_report.json containing final_residual. The
#               float-side clause (final_residual <= tol*(1+5e-4)) is
#               asserted by the pcg_e2e_verify gtest entry which reads the
#               report written here (CMake script math is integer-only).
#   rh_jacobi — run the jacobi e2e with --residual-history (np1, 32x32),
#               saving stdout to ${dir}/run_stdout.txt so the
#               ResidualHistoryE2ETest.JacobiFirstRowAndConvergedLine
#               verifier can compare the CSV last row against the
#               "Converged at iteration" log line at 6-digit print
#               precision.
#   badpath   — E6 (design §6): two variants of an unwritable
#               --residual-history path (nonexistent parent directory, and
#               a path whose parent is a regular file). Both must exit
#               non-zero with an [ERROR] log line (HYPOS_ERROR macro —
#               the logger prints "[ERROR]").
#   help      — W2 (design §6): `--help` must list the pcg solver value
#               (worded as pipelined, explicitly not preconditioned CG)
#               and the --residual-history option.
#
# Any violation terminates with FATAL_ERROR (non-zero exit).

if(NOT DEFINED exe OR NOT DEFINED dir OR NOT DEFINED mode)
    message(FATAL_ERROR
        "usage: cmake -Dexe=<hypos> -Ddir=<dir> -Dmode=<pcg_e2e|rh_jacobi|badpath|help> "
        "-P VerifyAr009Cli.cmake")
endif()

if(mode STREQUAL "pcg_e2e")
    file(MAKE_DIRECTORY "${dir}")
    execute_process(
        COMMAND "${exe}" --solver pcg --nx 32 --ny 32 --tol 1e-6
                --output-dir "${dir}"
        WORKING_DIRECTORY "${dir}"
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "ar009 cli verify (pcg_e2e): hypos exited ${rc}; output: ${out}; error: ${err}")
    endif()
    if(NOT out MATCHES "PipelinedCG finished")
        message(FATAL_ERROR "ar009 cli verify (pcg_e2e): output lacks the \"PipelinedCG finished\" line; output: ${out}")
    endif()
    if(NOT EXISTS "${dir}/performance_report.json")
        message(FATAL_ERROR "ar009 cli verify (pcg_e2e): performance_report.json missing in ${dir}; output: ${out}; error: ${err}")
    endif()
    file(READ "${dir}/performance_report.json" json)
    if(NOT json MATCHES "\"final_residual\":")
        message(FATAL_ERROR "ar009 cli verify (pcg_e2e): report JSON lacks final_residual; json: ${json}")
    endif()

elseif(mode STREQUAL "rh_jacobi")
    file(MAKE_DIRECTORY "${dir}")
    execute_process(
        COMMAND "${exe}" --solver jacobi --nx 32 --ny 32 --tol 1e-6
                --output-dir "${dir}"
                --residual-history "${dir}/residual_jacobi.csv"
        WORKING_DIRECTORY "${dir}"
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "ar009 cli verify (rh_jacobi): hypos exited ${rc}; output: ${out}; error: ${err}")
    endif()
    if(NOT EXISTS "${dir}/residual_jacobi.csv")
        message(FATAL_ERROR "ar009 cli verify (rh_jacobi): residual_jacobi.csv missing in ${dir}; output: ${out}; error: ${err}")
    endif()
    if(NOT EXISTS "${dir}/performance_report.json")
        message(FATAL_ERROR "ar009 cli verify (rh_jacobi): performance_report.json missing in ${dir}; output: ${out}; error: ${err}")
    endif()
    # Preserve stdout for the gtest-side "Converged at iteration" line
    # comparison (E2①, 6-digit print precision).
    file(WRITE "${dir}/run_stdout.txt" "${out}")

elseif(mode STREQUAL "badpath")
    file(MAKE_DIRECTORY "${dir}")
    # Variant 1: parent directory does not exist.
    execute_process(
        COMMAND "${exe}" --nx 16 --ny 16 --max-iter 2
                --output-dir "${dir}"
                --residual-history "${dir}/no_such_dir/f.csv"
        WORKING_DIRECTORY "${dir}"
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(rc EQUAL 0)
        message(FATAL_ERROR "ar009 cli verify (badpath): nonexistent-dir variant exited 0 (must fail); output: ${out}; error: ${err}")
    endif()
    if(NOT out MATCHES "\\[ERROR\\]")
        message(FATAL_ERROR "ar009 cli verify (badpath): nonexistent-dir variant lacks an [ERROR] line; output: ${out}; error: ${err}")
    endif()
    # Variant 2: parent path is a regular file (illegal path).
    file(WRITE "${dir}/anchor_file.txt" "x")
    execute_process(
        COMMAND "${exe}" --nx 16 --ny 16 --max-iter 2
                --output-dir "${dir}"
                --residual-history "${dir}/anchor_file.txt/f.csv"
        WORKING_DIRECTORY "${dir}"
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(rc EQUAL 0)
        message(FATAL_ERROR "ar009 cli verify (badpath): file-as-parent variant exited 0 (must fail); output: ${out}; error: ${err}")
    endif()
    if(NOT out MATCHES "\\[ERROR\\]")
        message(FATAL_ERROR "ar009 cli verify (badpath): file-as-parent variant lacks an [ERROR] line; output: ${out}; error: ${err}")
    endif()

elseif(mode STREQUAL "help")
    execute_process(
        COMMAND "${exe}" --help
        RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "ar009 cli verify (help): hypos --help exited ${rc}; output: ${out}; error: ${err}")
    endif()
    if(NOT out MATCHES "pcg")
        message(FATAL_ERROR "ar009 cli verify (help): --help solver line lacks pcg; output: ${out}")
    endif()
    if(NOT out MATCHES "pipelined")
        message(FATAL_ERROR "ar009 cli verify (help): --help pcg wording lacks \"pipelined\"; output: ${out}")
    endif()
    if(NOT out MATCHES "residual-history")
        message(FATAL_ERROR "ar009 cli verify (help): --help output lacks --residual-history; output: ${out}")
    endif()

else()
    message(FATAL_ERROR "ar009 cli verify: unknown mode '${mode}' (expected pcg_e2e|rh_jacobi|badpath|help)")
endif()

message(STATUS "ar009 cli verify (${mode}): ok")
