# VerifyAr005Output.cmake — AR005 B1 (design §6.3): assert that a np=1 e2e
# run with --output-format vtk|binary produced exactly the expected file
# set, and verify the on-disk payload contracts beyond mere existence:
#
#   vtk:    the per-rank piece must be appended raw XML — it contains
#           header_type="UInt64" and format="appended" (no format="ascii"),
#           and the little-endian UInt64 length header right after the
#           AppendedData '_' mark equals nx*ny*nz*8 bytes, with the file
#           ending within 64 bytes after the payload (closing tags only).
#   binary: the per-rank shard size is exactly 56 + nx*ny*nz*8 bytes, and
#           the first three 8-byte little-endian header fields decode to
#           nx, ny, nz.
#
# Usage:
#   cmake -Ddir=<dir> -Dexpected="f1;f2" -Dfmt=vtk|binary -Dpiece=<file>
#         -Dnx=16 -Dny=16 -Dnz=1 -P VerifyAr005Output.cmake
#
# Fails (non-zero exit) on any mismatch.

if(NOT DEFINED dir OR NOT DEFINED expected OR NOT DEFINED fmt
        OR NOT DEFINED piece OR NOT DEFINED nx OR NOT DEFINED ny OR NOT DEFINED nz)
    message(FATAL_ERROR
        "usage: cmake -Ddir=<dir> -Dexpected=<f1;f2> -Dfmt=vtk|binary -Dpiece=<file> "
        "-Dnx=.. -Dny=.. -Dnz=.. -P VerifyAr005Output.cmake")
endif()

if(NOT EXISTS "${dir}")
    message(FATAL_ERROR "output dir does not exist: ${dir}")
endif()

# --- exact file-set check (mirrors VerifyMpiioOutput.cmake) ---
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
        "unexpected ${fmt} output files in ${dir}\n"
        "  actual:   ${actual_names}\n"
        "  expected: ${expected_names}")
endif()

set(piece_path "${dir}/${piece}")
if(NOT EXISTS "${piece_path}")
    message(FATAL_ERROR "piece file does not exist: ${piece_path}")
endif()

file(SIZE "${piece_path}" sz)
math(EXPR data_bytes "${nx} * ${ny} * ${nz} * 8")

if(fmt STREQUAL "vtk")
    # The XML prologue is plain ASCII and contains no NUL bytes, so a text
    # read yields the prologue verbatim; it stops at the first NUL of the
    # binary payload, which is fine for the structural checks below.
    file(READ "${piece_path}" prologue)
    foreach(needle "header_type=\"UInt64\"" "format=\"appended\""
            "<AppendedData encoding=\"raw\">")
        if(NOT prologue MATCHES "${needle}")
            message(FATAL_ERROR "piece ${piece} lacks '${needle}' (not appended raw?)")
        endif()
    endforeach()
    if(prologue MATCHES "format=\"ascii\"")
        message(FATAL_ERROR "piece ${piece} still declares format=\"ascii\"")
    endif()
    # The '_' mark follows the AppendedData open tag. Note the XML prologue
    # contains other underscores (header_type!), so search relative to the
    # AppendedData tag instead of taking the first '_' in the file.
    set(ad_tag "<AppendedData encoding=\"raw\">")
    string(FIND "${prologue}" "${ad_tag}" ad_idx)
    if(ad_idx LESS 0)
        message(FATAL_ERROR "piece ${piece}: AppendedData open tag not found")
    endif()
    string(LENGTH "${ad_tag}" ad_len)
    math(EXPR after "${ad_idx} + ${ad_len}")
    string(SUBSTRING "${prologue}" "${after}" -1 rest)
    string(FIND "${rest}" "_" rel_idx)
    if(rel_idx LESS 0)
        message(FATAL_ERROR "piece ${piece}: AppendedData '_' mark not found")
    endif()
    math(EXPR mark_idx "${after} + ${rel_idx}")

    file(READ "${piece_path}" hex HEX)
    # The HEX conversion inserts a newline every 32 hex chars; strip them so
    # string offsets map 1:2 onto byte offsets.
    string(REGEX REPLACE "\n" "" hex "${hex}")
    math(EXPR hdr_off "${mark_idx} * 2 + 2")
    string(SUBSTRING "${hex}" "${hdr_off}" 16 len_le)
    string(SUBSTRING "${len_le}" 8 8 len_hi)
    if(NOT len_hi STREQUAL "00000000")
        message(FATAL_ERROR "piece ${piece}: payload length >= 2^32 (unsupported)")
    endif()
    string(SUBSTRING "${len_le}" 0 2 b0)
    string(SUBSTRING "${len_le}" 2 2 b1)
    string(SUBSTRING "${len_le}" 4 2 b2)
    string(SUBSTRING "${len_le}" 6 2 b3)
    math(EXPR len_hdr "0x${b3}${b2}${b1}${b0}")
    if(NOT len_hdr EQUAL data_bytes)
        message(FATAL_ERROR
            "piece ${piece}: appended length header ${len_hdr} != nx*ny*nz*8 (${data_bytes})")
    endif()
    # Payload must end within 64 bytes of the file end (closing tags only).
    math(EXPR tail "${sz} - (${mark_idx} + 9 + ${len_hdr})")
    if(tail LESS 1 OR tail GREATER 63)
        message(FATAL_ERROR
            "piece ${piece}: ${tail} bytes after appended payload (expected closing tags only, <64)")
    endif()
elseif(fmt STREQUAL "binary")
    math(EXPR want "56 + ${data_bytes}")
    if(NOT sz EQUAL want)
        message(FATAL_ERROR "shard ${piece}: size ${sz} != 56 + nx*ny*nz*8 (${want})")
    endif()
    file(READ "${piece_path}" hex HEX)
    string(REGEX REPLACE "\n" "" hex "${hex}")
    # Header fields: nx, ny, nz at byte offsets 0, 8, 16 (little-endian).
    set(fields nx;ny;nz)
    set(values ${nx};${ny};${nz})
    set(offs 0;8;16)
    set(idx 0)
    foreach(fname IN LISTS fields)
        list(GET values ${idx} fval)
        list(GET offs ${idx} foff)
        math(EXPR hoff "${foff} * 2")
        string(SUBSTRING "${hex}" "${hoff}" 16 f_le)
        string(SUBSTRING "${f_le}" 0 2 b0)
        string(SUBSTRING "${f_le}" 2 2 b1)
        string(SUBSTRING "${f_le}" 4 2 b2)
        string(SUBSTRING "${f_le}" 6 2 b3)
        string(SUBSTRING "${f_le}" 8 8 f_hi)
        if(NOT f_hi STREQUAL "00000000")
            message(FATAL_ERROR "shard ${piece}: header ${fname} >= 2^32 (unsupported)")
        endif()
        math(EXPR got "0x${b3}${b2}${b1}${b0}")
        if(NOT got EQUAL fval)
            message(FATAL_ERROR "shard ${piece}: header ${fname}=${got}, expected ${fval}")
        endif()
        math(EXPR idx "${idx} + 1")
    endforeach()
else()
    message(FATAL_ERROR "fmt must be vtk or binary, got: ${fmt}")
endif()

message(STATUS "${fmt} output verified (${dir}): ${actual_names}")
