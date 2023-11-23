#
# Nvtx can be used even when not using GPU, for a pure CPU run.
# All we need to perfrom tracing is to run on a platform where cuda driver
# and nsys profile/tracing tools are installed
#

#
# Nvtx is already shipped with nvcc toolkit and nvhpc compiler, but we want to be able to
# instrument code event for pure CPU build and run; so we need a way to provide nvtx when
# not building with nvcc or nvc++
#

# 1. if using an Nvidia compiler, we don't need to download nvtx sources (they are already available)
# 2. else we need to download nvtx sources to enable NVTX annotations

option (EULER2D_NVTX_ANNOTATION_ENABLED "build with nvtx annotations enabled (default OFF)" OFF)

#
# Option to use git (instead of tarball release) for downloading nvtx
#
option(EULER2D_NVTX_USE_GIT "Turn ON if you want to use git (instead of archive file) to download nvtx sources (default: ON)" ON)

set(EULER2D_USE_NVTX_FROM_GITHUB FALSE)
set(EULER2D_USE_NVTX3 OFF)

# If use requested NVTX annotations, we need to check which compiler is used.
# If compiler is not an nvidia compiler, we need to download nvtx sources
if(EULER2D_NVTX_ANNOTATION_ENABLED)

  # check if compiler is nvcc or nvc++, nvtx is available
  if (CMAKE_CXX_COMPILER_ID MATCHES "NVHPC" OR KOKKOS_CXX_COMPILER_ID MATCHES "NVIDIA")
    message("[euler2d / nvtx] NVTX found, using a nvidia compiler")
    set(EULER2D_NVTX3_FOUND TRUE)
  else()

    message("[euler2d / nvtx] Downloading nvtx sources")

    set_property(DIRECTORY PROPERTY EP_BASE ${CMAKE_BINARY_DIR}/external)

    include (FetchContent)

    if (EULER2D_NVTX_USE_GIT)
      FetchContent_Declare( nvtx_external
        GIT_REPOSITORY https://github.com/nvidia/nvtx.git
        GIT_TAG release-v3
        SOURCE_SUBDIR c
        )
    else()
      FetchContent_Declare( nvtx_external
        URL https://github.com/NVIDIA/NVTX/archive/refs/tags/v3.1.0.tar.gz
        )
    endif()

    # Import nvtx targets (download, and call add_subdirectory)
    FetchContent_MakeAvailable(nvtx_external)

    if(TARGET nvtx3-c)
      message("[euler2d / nvtx] NVTX found (using FetchContent)")
      set(EULER2D_NVTX3_FOUND True)
    else()
      message("[euler2d / nvtx] we shouldn't be here. We've just integrated nvtx build into euler2d build !")
    endif()

    set(EULER2D_USE_NVTX_FROM_GITHUB TRUE)

  endif()

endif(EULER2D_NVTX_ANNOTATION_ENABLED)

if(EULER2D_NVTX3_FOUND)
  set(EULER2D_USE_NVTX3 ON)
endif()

if (EULER2D_NVTX_ANNOTATION_ENABLED AND EULER2D_NVTX3_FOUND)
  add_compile_options(-DEULER2D_NVTX_ANNOTATION_ENABLED)
endif(EULER2D_NVTX_ANNOTATION_ENABLED AND EULER2D_NVTX3_FOUND)
