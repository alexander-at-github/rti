cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

# Set compiler
#set(ENV{CXX} "/usr/bin/clang++-6.0")
#set(ENV{CXX} "/usr/bin/clang")
# REMARK: setting ENV{CXX} only works when the project function/macro is called afterwards!
## TODO: Which setting overrides the other? ENV or CACHE?

project(
  rti-superbuild
  VERSION 0.1
  LANGUAGES CXX C
  )

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

##########################################################################################
### If one wants to uses an Intel SPMD Compiler, the one has to set that here
##########################################################################################
set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "Set whether or not you use Intel's SPMD Compiler")
##########################################################################################

set(M_DEPENDENCIES_DIR ${CMAKE_SOURCE_DIR}/dependencies)
set(STAGED_INSTALL_PREFIX ${M_DEPENDENCIES_DIR}/stage) # will be used in ./external/upstream/
message(STATUS "${PROJECT_NAME} external project's staged install path: ${STAGED_INSTALL_PREFIX}")

list(APPEND BOOST_COMPONENTS_REQUIRED "log" "system" "timer" "iostreams" "filesystem")

include(ExternalProject)
add_subdirectory(external/upstream)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "NOT CMAKE_BUILD_TYPE holds")
  #set(CMAKE_BUILD_TYPE Debug)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Set build type to \"Debug\" or \"Release\"" FORCE)
endif()
message(STATUS "CMAKE_BUILD_TYPE == ${CMAKE_BUILD_TYPE}")


# -Wno-unused-parameter
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto") #"-DNDEBUG -O3 -march=native -flto")
#set(CMAKE_CXX_FLAGS_RELEASE "-O")
#set(CMAKE_CXX_FLAGS_RELEASE "-pg -O") # arguments for gprof
#set(CMAKE_CXX_FLAGS_RELEASE "-O -fprofile-arcs -ftest-coverage") # arguments for gcov

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}") # CACHE STRING "cmake c++ flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "cmake c++ flags" FORCE)
  message(STATUS "CMAKE_BUILD_TYPE == Debug and CMAKE_CXX_FLAGS == ${CMAKE_CXX_FLAGS}")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}") # CACHE STRING "cmake c++ flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "cmake c++ flags" FORCE)
  message(STATUS "CMAKE_BUILD_TYPE == Release and CMAKE_CXX_FLAGS == ${CMAKE_CXX_FLAGS}")
endif()
message(STATUS "CMAKE_CXX_FLAGS == ${CMAKE_CXX_FLAGS}")

externalproject_add(rti
  DEPENDS
		#embree-external
    #vtk-external
    #
    # optional:
    #gtest-external
    #abseil-external
    #
    # The project depends on embree-external and other external libraries. We do not add it here,
    # because otherwise the build system is very slow when rebuilding the project.
    # (The build configuration of these dependecies redo a lot of steps). You are adviced
    # to build these targets once yourself.
  SOURCE_DIR
    ${CMAKE_CURRENT_LIST_DIR}/src
  # Eventhough this is a superbuild, we want the final binary in the build directory
  INSTALL_DIR
    ${CMAKE_BINARY_DIR}
  #LIST_SEPARATOR ;
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
    -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
    #
    -DBOOST_ROOT:PATH=${BOOST_ROOT}
    -DGMSH_DIR=${GMSH_DIR}
    -DEMBREE_DIR=${EMBREE_DIR}
    -DVTK_DIR=${VTK_DIR}
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
    -DCMAKE_INCLUDE_PATH:PATH=${CMAKE_INCLUDE_PATH}
    -DCMAKE_LIBRARY_PATH:PATH=${CMAKE_LIBRARY_PATH}
    # To pass a list into externalproject_add() would need some yoga.
    # One can see an example here
    # https://stackoverflow.com/questions/45414507/pass-a-list-of-prefix-paths-to-externalproject-add-in-cmake-args
    #-D BOOST_COMPONENTS_REQUIRED="${BOOST_COMPONENTS_REQUIRED}"
  BUILD_ALWAYS
    1
  #INSTALL_COMMAND
  #  ""
)

externalproject_get_property(rti SOURCE_DIR)
set(RTI_SRC_DIR ${SOURCE_DIR})
unset(SOURCE_DIR)

externalproject_add(rti-test
  DEPENDS
  #gtest-external
  #boost-external
  #gnuplot-iostream-external
  SOURCE_DIR
  ${CMAKE_CURRENT_LIST_DIR}/test
  INSTALL_DIR
  ${CMAKE_BINARY_DIR}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
    -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
    #
    -DBOOST_ROOT:PATH=${BOOST_ROOT}
    #-DGMSH_DIR=${GMSH_DIR}
    #-DEMBREE_DIR=${EMBREE_DIR}
    #-DVTK_DIR=${VTK_DIR}
    #
    -DRTI_SRC_DIR=${RTI_SRC_DIR}
    #
    -DGTEST_SRC_DIR=${GTEST_SRC_DIR}
    -DGTEST_BIN_DIR=${GTEST_BIN_DIR}
    -DGPLOT_IO_DIR=${GPLOT_IO_DIR}
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
    -DCMAKE_INCLUDE_PATH:PATH=${CMAKE_INCLUDE_PATH}
    -DCMAKE_LIBRARY_PATH:PATH=${CMAKE_LIBRARY_PATH}
    #-D BOOST_COMPONENTS_REQUIRED="${BOOST_COMPONENTS_REQUIRED}"
  BUILD_ALWAYS
    1
  )

externalproject_add(rti-bench
  DEPENDS
  #gbench-external
  SOURCE_DIR
  ${CMAKE_CURRENT_LIST_DIR}/bench
  INSTALL_DIR
  ${CMAKE_BINARY_DIR}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
    -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
    #
    -DBOOST_ROOT:PATH=${BOOST_ROOT}
    #-DGMSH_DIR=${GMSH_DIR}
    -DEMBREE_DIR=${EMBREE_DIR}
    #-DVTK_DIR=${VTK_DIR}
    #
    -DRTI_SRC_DIR=${RTI_SRC_DIR}
    #
    -DGBENCH_SRC_DIR=${GBENCH_SRC_DIR}
    -DGBENCH_BIN_DIR=${GBENCH_BIN_DIR}
    -DGBENCH_INST_DIR=${GBENCH_INST_DIR}
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
    -DCMAKE_INCLUDE_PATH:PATH=${CMAKE_INCLUDE_PATH}
    -DCMAKE_LIBRARY_PATH:PATH=${CMAKE_LIBRARY_PATH}
    #-D BOOST_COMPONENTS_REQUIRED="${BOOST_COMPONENTS_REQUIRED}"
  BUILD_ALWAYS
    1
  )

externalproject_add(rti-experimental
  DEPENDS
  #
  SOURCE_DIR
  ${CMAKE_CURRENT_LIST_DIR}/experimental
  INSTALL_DIR
  ${CMAKE_BINARY_DIR}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
    -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
    #
    -DBOOST_ROOT:PATH=${BOOST_ROOT}
    -DEMBREE_DIR=${EMBREE_DIR}
    #-DVTK_DIR=${VTK_DIR}
    #
    -DRTI_SRC_DIR=${RTI_SRC_DIR}
  CMAKE_CACHE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
    -DCMAKE_INCLUDE_PATH:PATH=${CMAKE_INCLUDE_PATH}
    -DCMAKE_LIBRARY_PATH:PATH=${CMAKE_LIBRARY_PATH}
  BUILD_ALWAYS
    1
  )
