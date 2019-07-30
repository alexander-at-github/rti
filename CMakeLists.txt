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

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

##########################################################################################
### If one wants to uses an Intel SPMD Compiler, the one has to set that here
##########################################################################################
set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "Set whether or not you use Intel's SPMD Compiler")
##########################################################################################

set(M_DEPENDENCIES_DIR ${CMAKE_SOURCE_DIR}/dependencies)
set(STAGED_INSTALL_PREFIX ${M_DEPENDENCIES_DIR}/stage) # will be used in external/upstream
message(STATUS "${PROJECT_NAME} external project's staged install path: ${STAGED_INSTALL_PREFIX}")

list(APPEND BOOST_COMPONENTS_REQUIRED "log" "system" "timer")

include(ExternalProject)
add_subdirectory(external/upstream)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
# -Wno-unused-parameter
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -march=native")

externalproject_add(rti
  DEPENDS
    #boost-external
    #gmsh-external
		#embree-external
    #vtk-external
    #
    # The project depends on gmsh-external and embree-external. We do not add it here,
    # because otherwise the build system is very slow when rebuilding the project.
    # (The build configuration of these dependecies redo a lot of steps). You are adviced
    # to build these two targets once yourself with the following commands:
    # $ cmake --build . --target embree-external
    # $ cmake --build . --target gmsh-external
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
