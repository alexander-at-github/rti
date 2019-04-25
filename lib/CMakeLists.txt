cmake_minimum_required(VERSION 3.7)
message(STATUS "CMAKE_ROOT=${CMAKE_ROOT}")

## Compiler Options
set(CMAKE_C_COMPILER               "/usr/bin/gcc" CACHE FILEPATH "C compiler set to gcc")
set(CMAKE_CXX_COMPILER             "/usr/bin/g++" CACHE FILEPATH "CXX compiler set to g++")
set(CMAKE_CXX_FLAGS                " -Wall -std=c++14 -DBOOST_LOG_DYN_LINK -DBOOST_SYSTEM_NO_DEPRECATED" CACHE STRING "compile flags" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG          " -O0 -g" CACHE STRING "Set CXX debug flags")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO " -O2 -g" CACHE STRING "Set some CXX flags")
set(CMAKE_CXX_FLAGS_RELEASE        " -O3 -DNDEBUG" CACHE STRING "Set release CXX flags")

## Set build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

## General Set-Up
#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#    set(HPTCAD_TOOLSET "clang")
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    set(HPTCAD_TOOLSET "gcc")
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#    set(HPTCAD_TOOLSET "intel-linux")
#endif()
#message(status "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
#message(status "CMAKE_CURRENT_LIST_DIR=${CMAKE_CURRENT_LIST_DIR}")
message("CMAKE_C_COMPILER=" ${CMAKE_C_COMPILER})
message("CMAKE_CXX_COMPILER=" ${CMAKE_CXX_COMPILER})
message("CMAKE_CXX_COMPILER_ID=" ${CMAKE_CXX_COMPILER_ID})
message("CMAKE_BUILD_TYPE=" ${CMAKE_BUILD_TYPE})
#message("HPTCAD_TOOLSET=" ${HPTCAD_TOOLSET})

## Threading Building Blocks
message("#############")
message("# Build TBB #")
message("#############")
set(TBB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/tbb")
include("${TBB_ROOT}/cmake/TBBBuild.cmake")
tbb_build(TBB_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/tbb CONFIG_DIR TBB_DIR)
find_package(TBB REQUIRED)
#message(STATUS "${TBB_IMPORTED_TARGETS}")
#target_link_libraries(<executable> TBB::tbb)
message("# Build TBB Done #")

## Optionally we can use out-of-OS GLFW (https://github.com/glfw/glfw/)

### Embree 3
message("#################")
message("# Build Embree3 #")
message("#################")
set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "We do not use the Intel SPMD Compiler")
#set(EMBREE_MAX_ISA AVX2) #SSE2, SSE4.2, AVX, AVX2, AVX512KNL, AVX512SKX
add_subdirectory(embree)
message("# Build Embree3 Done #")
