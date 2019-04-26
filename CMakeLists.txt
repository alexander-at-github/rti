cmake_minimum_required(VERSION 3.7)
project(rti)


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


## Add libraries
#
# Compile Embree as a static library.
# If we want to use Embree as dynamic library, then one needs to comment
# the following line out and add the location of the dynamic library
# (<build-folder>/lib/embree) to the environment variable LD_LIBRARY_PATH.
#set(EMBREE_STATIC_LIB ON CACHE BOOL "Compile Embree static library")
add_subdirectory(lib)
# Add package Embree 3.5
set(embree_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/embree")
find_package(embree 3.5 REQUIRED)
## Add Embree headers to include directories
include_directories("${EMBREE_INCLUDE_DIRS}")

## Boost
find_package(Boost COMPONENTS log REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

file(GLOB SOURCES "src/**/*.cpp")
add_executable(TestRTIntersect ${SOURCES})
#message(STATUS "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
# Add embree as dependecy to our executable.
add_dependencies(TestRTIntersect embree)
#message(STATUS "EMBREE_LIBRARIES=${EMBREE_LIBRARIES}")
#message(STATUS "EMBREE_LIBRARY=${EMBREE_LIBRARY}")
# [Question] Exactly which file do we link here?
target_link_libraries(TestRTIntersect embree)
#target_link_libraries(TestRTIntersect ${EMBREE_LIBRARY})
target_link_libraries(TestRTIntersect "${CMAKE_CURRENT_SOURCE_DIR}/lib/gmsh/install/lib/libgmsh.a")
## One needs to link lapack and blas (in that order) for gmesh!
target_link_libraries(TestRTIntersect lapack blas)
set(OCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/oce/install/oce/lib/oce-0.18")
find_package(OCE REQUIRED)
target_link_libraries(TestRTIntersect ${OCE_LIBRARIES})
target_link_libraries(TestRTIntersect ${Boost_LOG_LIBRARY})


## gmsh
#set(gmsh_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gmsh/build/install")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/lib/gmsh/install/include")
