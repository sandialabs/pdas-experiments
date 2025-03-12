# Largely copied from pressiodemoapps, to make sure we're compatible

cmake_minimum_required(VERSION 3.18.0)
project(pdas-experiments CXX)

# fetch pressio-ops
if(PRESSIO_OPS_INCLUDE_DIR)
  message("-- Using ${PRESSIO_OPS_INCLUDE_DIR}")
  include_directories(SYSTEM ${PRESSIO_OPS_INCLUDE_DIR})

else()

  include(FetchContent)
  message("-- PRESSIO_OPS_INCLUDE_DIR is not defined, so fetching pressio-ops...")

  FetchContent_Declare(
    pressio-ops
    GIT_REPOSITORY https://github.com/Pressio/pressio-ops.git
    GIT_TAG        develop
    GIT_PROGRESS   TRUE
    GIT_SHALLOW    TRUE
  )

  FetchContent_MakeAvailable(pressio-ops)
  set(PRESSIO_OPS_INCLUDE_DIR ${pressio-ops_SOURCE_DIR}/include)
endif()

# fetch pressio-log
if(PRESSIO_LOG_INCLUDE_DIR)
  message("-- Using ${PRESSIO_LOG_INCLUDE_DIR}")
  include_directories(SYSTEM ${PRESSIO_LOG_INCLUDE_DIR})

else()

  include(FetchContent)
  message("-- PRESSIO_LOG_INCLUDE_DIR is not defined, so fetching pressio-log...")

  FetchContent_Declare(
    pressio-log
    GIT_REPOSITORY https://github.com/Pressio/pressio-log.git
    GIT_TAG        main
    GIT_PROGRESS   TRUE
    GIT_SHALLOW    TRUE
  )

  FetchContent_MakeAvailable(pressio-log)
  set(PRESSIO_LOG_INCLUDE_DIR ${pressio-log_SOURCE_DIR}/include)
endif()

#=====================================================================
# we need c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORT_CPP17)
if(NOT COMPILER_SUPPORT_CPP17)
  message(FATAL_ERROR "Compiler does not support -std=c++17. This is required.")
endif()
set(PRESSIO_ENABLE_CXX17 ON)

#=====================================================================
# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed.
    Please make a new directory (e.g. build directory) and run CMake from there.")
endif()

# default to release if build type is empty
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "${CMAKE_BUILD_TYPE} is not specified, defaulting to Release.")
  set(CMAKE_BUILD_TYPE "Release")
endif()

# convert cmake build type to lower string
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if( NOT cmake_build_type_tolower STREQUAL "debug"
    AND NOT cmake_build_type_tolower STREQUAL "release")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
    Allowed values are Debug, Release (case-insensitive).")
endif()

#=====================================================================
# We're not going to mess with OpenMP or Python bindings for now
# Just check that pressio, pressiodemoapps, and pressio-schwarz are specified correctly

if(NOT PDA_SOURCE)
  message(FATAL_ERROR "Must specify pressiodemoapps source directory via PDA_SOURCE")
else()
  if(NOT EXISTS ${PDA_SOURCE})
    message(NOTICE "PDA_SOURCE: ${PDA_SOURCE}")
    message(FATAL_ERROR "PDA_SOURCE does not exist or is not an absolute path")
  endif()
endif()

if(NOT PRESSIO_SOURCE)
  message(FATAL_ERROR "Must specify pressio source directory via PRESSIO_SOURCE")
else()
  if(NOT EXISTS ${PRESSIO_SOURCE})
    message(NOTICE "PRESSIO_SOURCE: ${PRESSIO_SOURCE}")
    message(FATAL_ERROR "PRESSIO_SOURCE does not exist or is not an absolute path")
  endif()
endif()

if(NOT PSCHWARZ_SOURCE)
  message(FATAL_ERROR "Must specify pressio-schwarz source directory via PSCHWARZ_SOURCE")
else()
  if(NOT EXISTS ${PSCHWARZ_SOURCE})
    message(NOTICE "PSCHWARZ_SOURCE: ${PSCHWARZ_SOURCE}")
    message(FATAL_ERROR "PSCHWARZ_SOURCE does not exist or is not an absolute path")
  endif()
endif()

# build YAML
add_subdirectory(
  ${CMAKE_CURRENT_SOURCE_DIR}/tpls/yaml-cpp-0.7.0
  ${CMAKE_CURRENT_BINARY_DIR}/yamlcpp)

# include demoapps, Schwarz
include_directories(
  ${PRESSIO_OPS_INCLUDE_DIR}
  ${PRESSIO_LOG_INCLUDE_DIR}
  ${PDA_SOURCE}/tpls/eigen3
  ${PDA_SOURCE}/include
  ${PRESSIO_SOURCE}/include
  ${PSCHWARZ_SOURCE}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/tpls/yaml-cpp-0.7.0/include)

add_compile_definitions(SCHWARZ_SAVE_TEMPDIR)
link_libraries(stdc++fs)
add_compile_definitions(PRESSIO_ENABLE_TPL_EIGEN PRESSIO_ENABLE_DEBUG_PRINT PRESSIO_ENABLE_INTERNAL_SPDLOG PRESSIO_ENABLE_LOGGING)

# serial executable
add_executable(runner_serial ${CMAKE_CURRENT_SOURCE_DIR}/src/runner.cc)
target_link_libraries(runner_serial PRIVATE yaml-cpp)

# parallel executable
add_executable(runner_omp ${CMAKE_CURRENT_SOURCE_DIR}/src/runner.cc)
target_link_libraries(runner_omp PRIVATE yaml-cpp)
find_package(OpenMP)
target_compile_definitions(runner_omp PRIVATE SCHWARZ_ENABLE_OMP)
target_link_libraries(runner_omp PRIVATE OpenMP::OpenMP_CXX pthread)
target_compile_options(runner_omp PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)

