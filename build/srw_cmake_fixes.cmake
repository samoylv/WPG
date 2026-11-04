# Fixes for the SRW CMake build, injected by build_srw.py with -DCMAKE_PROJECT_INCLUDE
# (runs right after project(SRW), SRW sources are not changed).
#
# 1. When FFTW is not found on the system, SRW builds it as an ExternalProject
#    (fftw_external for FFTW2 / OpenMP, fftw3_external for FFTW3), but the SRW
#    libraries do not depend on it. In a parallel build srwlpy is then linked before
#    libfftw.a exists ("No rule to make target stage/lib/libfftw.a").
# 2. C++ sources of the CUDA build include cufft.h, but only the CUDA sources get
#    the CUDA toolkit include directory.

if(USE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    include_directories(${CUDAToolkit_INCLUDE_DIRS})
endif()

function(_wpg_add_fftw_dependencies)
    foreach(ext fftw_external fftw3_external)
        if(TARGET ${ext})
            foreach(tgt srw srwlpy srwlclient)
                if(TARGET ${tgt})
                    add_dependencies(${tgt} ${ext})
                endif()
            endforeach()
        endif()
    endforeach()
endfunction()

# at the end of the top-level CMakeLists.txt, when all targets exist
cmake_language(DEFER DIRECTORY ${CMAKE_SOURCE_DIR} CALL _wpg_add_fftw_dependencies)
