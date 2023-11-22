cmake_minimum_required(VERSION 3.23)

#
# Allow find_package() to use <PackageName>_ROOT variables,
# either cmake variable or environment variables
# see https://cmake.org/cmake/help/latest/policy/CMP0074.html
#
if(NOT CMAKE_VERSION VERSION_LESS 3.12)
  cmake_policy(SET CMP0074 NEW)
endif()

project(abinit
  LANGUAGES C CXX Fortran
  HOMEPAGE_URL https://github.com/pkestene/euler2d_cudaFortran
  DESCRIPTION "Simple example of standard parallelism in Fortran.")

#
# Export compile command to json (for editors like emacs, clion, vim, etc...).
# It allows nice code editing features provided by LSP (Language Server Protocol)
#
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

#
# Prevent build from source directory
#
# make sure the user doesn't play dirty with symlinks
get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)

# disallow in-source builds
if("${srcdir}" STREQUAL "${bindir}")
  message("######################################################")
  message("# ${PROJECT_NAME} should not be configured and built in the source directory")
  message("# You must run cmake from a build directory.")
  message("# For example:")
  message("# mkdir _build ; cd _build")
  message("# run cmake from the build directory.")
  message("######################################################")
  message(FATAL_ERROR "Quitting configuration")
endif()

#
# Set default compile optimization flag
#
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(DEFAULT_BUILD_TYPE "Release")
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

include(cmake/compiler_setup.cmake)

add_executable(euler2d_stdpar_gpu
  m_nvtx.f90
  m_precision.f90
  m_constants.f90
  m_parameters.f90
  m_monitoring.f90
  m_utils.f90
  m_run.f90
  main.f90)
if(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC" OR CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  target_compile_options(euler2d_stdpar_gpu PRIVATE -stdpar=gpu -acc=gpu -gpu=ccnative -gpu=rdc,managed -Minform=warn -Minfo -Mpreprocess)
  target_link_options(euler2d_stdpar_gpu PRIVATE -stdpar=gpu -acc=gpu -gpu=rdc,managed)

  set(NVTX_LDFLAGS -L$ENV{NVHPC_ROOT}/cuda/lib64 -lnvToolsExt)
  target_link_options(euler2d_stdpar_gpu PRIVATE ${NVTX_LDFLAGS})

endif()
