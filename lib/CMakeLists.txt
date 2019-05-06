cmake_minimum_required(VERSION 3.7)

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

  ##  ### OCE
  ##  message("# Build OpenCascade Community Edition #")
  ##
  ##  file(MAKE_DIRECTORY "${OCE_INSTALL_PREFIX}")
  ##  set(OCE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/oce/install" CACHE STRING "Set OCE install directory")
  ##  ExternalProject_Add(OCE)
  ##  #ExternalProject_Add(OCE)
  ##  # cmake -DOCE_INSTALL_PREFIX=$(pwd)"/../install/oce" ..
  ##  # make -j8
  ##  # make install/strip
  ##
  ##  ### gmsh
  ##  # (export CASROOT=../oce/install/oce && echo $CASROOT && cmake -DCMAKE_INSTALL_PREFIX=$(pwd)/../install" ..)
  ##  # make
  ##  # make install
  ##  #
  ##  # additional information: https://gitlab.onelab.info/gmsh/gmsh/blob/master/README.txt

  ##  ##  (export CASROOT=../oce/install/oce && echo $CASROOT && cmake -DCMAKE_INSTALL_PREFIX="$(pwd)/../install" .. && make -j8 && make install) && rm -r /home/alexander/cdr/rti/build/* && (export CASROOT=../oce/install/oce && echo $CASROOT && cmake -DCMAKE_INSTALL_PREFIX="$(pwd)/../install" -DENABLE_FLTK=0 -DDEFAULT=0 -DENABLE_BUILD_LIB=1 -DENABLE_POST=1 -DENABLE_PARSER=1 .. && make lib -j8 && make install/fast)

  ##  ## (export CASROOT=../oce/install/oce && echo $CASROOT && cmake -DCMAKE_INSTALL_PREFIX="$(pwd)/../install" -DENABLE_FLTK=0 -DCMAKE_BUILD_TYPE=Debug -DDEFAULT=0 -DENABLE_BUILD_LIB=1 -DENABLE_BUILD_SHARED=1 -DENABLE_POST=1 -DENABLE_PARSER=1 ..)
