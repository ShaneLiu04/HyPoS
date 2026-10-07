# Compiler-specific warning, optimization and sanitizer options for HyPoS.
# Included from the top-level CMakeLists.txt after project().

if(MSVC)
    message(FATAL_ERROR
        "HyPoS targets POSIX/Linux (MPI + posix_memalign + OpenMP); "
        "MSVC is not a supported compiler.")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
    add_compile_options(
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-DNDEBUG>
        $<$<CONFIG:Debug>:-O0>
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Debug>:-fsanitize=address,undefined>
    )
    add_link_options($<$<CONFIG:Debug>:-fsanitize=address,undefined>)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(-Wall -Wextra)
    add_compile_options(
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-xHost>
        $<$<CONFIG:Release>:-DNDEBUG>
        $<$<CONFIG:Debug>:-O0>
        $<$<CONFIG:Debug>:-g>
    )
else()
    message(FATAL_ERROR
        "Unsupported compiler '${CMAKE_CXX_COMPILER_ID}': "
        "HyPoS supports GNU, Clang and Intel.")
endif()
