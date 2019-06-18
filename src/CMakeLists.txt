cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(
  rti
  VERSION 0.1
  LANGUAGES CXX
  )

message(STATUS "CMAKE_CXX_FLAGS == ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_BUILD_TYPE == '${CMAKE_BUILD_TYPE}'")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "CMAKE_CXX_FLAGS_DEBUG == ${CMAKE_CXX_FLAGS_DEBUG}")
endif()
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "CMAKE_CXX_FLAGS_RELEASE == ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

set(BUILD_SHARED_LIBS ON) # A particularity of Boost 1.70

find_package(Boost REQUIRED COMPONENTS log system timer # "${BOOST_COMPONENTS_REQUIRED}"
  PATHS ${BOOST_ROOT}
  NO_DEFAULT_PATH
  )
find_package(embree 3 REQUIRED
  PATHS ${EMBREE_DIR}
  NO_DEFAULT_PATH
  )
# Gmsh does not provide a config of a module for find_package().
# We find the files ourselfs.
find_path(GMSH_INCLUDE_DIR NAMES gmsh.h
  PATHS "${GMSH_DIR}/include"
	NO_DEFAULT_PATH
	)
find_library(GMSH_LIBRARY NAMES gmsh libgmsh
  PATHS "${GMSH_DIR}/lib"
  NO_DEFAULT_PATH
)
# TBB
message(STATUS "TBB_DIR == ${TBB_DIR}")
find_path(TBB_INCLUDE_DIR NAMES tbb/task_scheduler_init.h
  PATHS "${TBB_DIR}/include" # ${TBB_DIR}/include/tbb/task_scheduler_init.h
  NO_DEFAULT_PATH
)
find_library(TBB_LIBRARY NAMES tbb libtbb
  PATHS "${TBB_DIR}/build/lib_release" # ${TBB_DIR}/build/lib_release/libtbb.so
  NO_DEFAULT_PATH
  )
find_library(TBB_MALLOC_LIBRARY NAMES tbbmalloc libtbbmalloc
  PATHS "${TBB_DIR}/build/lib_release"
  NO_DEFAULT_PATH
  )
message(STATUS "TBB_INCLUDE_DIR == ${TBB_INCLUDE_DIR}")
message(STATUS "TBB_LIBRARY == ${TBB_LIBRARY}")
message(STATUS "TBB_MALLOC_LIBRARY == ${TBB_MALLOC_LIBRARY}")
message(STATUS "TODO: use tbb-malloc")


# One needs to set the RPATH, otherwise the installed executable does not have dynamic
# link libraries set correctly.
# CMAKE_INSTALL_RPATH_USE_LINK_PATH needs to be set before the executables are added.
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# We set this property later on the executable itself.


# Executable for testing the ray tracing code
add_executable(rtiex rti/main_testing_rt.cpp)
target_include_directories(rtiex
  PRIVATE
	  ${CMAKE_CURRENT_LIST_DIR}
    ${EMBREE_INCLUDE_DIRS}
	  ${GMSH_INCLUDE_DIR}
    ${TBB_INCLUDE_DIR}
  )
target_link_libraries(rtiex
  PRIVATE
    Boost::log
    Boost::system
    Boost::timer
    ${EMBREE_LIBRARY}
		${GMSH_LIBRARY}
    ${TBB_LIBRARY}
  )
# Copy executable straight into the directory given by the variable CMAKE_INSTALL_PREFIX.
install(TARGETS rtiex
  DESTINATION "."
  # If one wants to install the executable somewhere else, then one should consider
  # to include GNUInstallDirs and the variable CMAKE_INSTALL_BINDIR.
  )
# One needs to set the RPATH, otherwise the installed executable does not have dynamic
# link libraries set correctly.
set_target_properties(rtiex
  PROPERTIES
  INSTALL_RPATH_USE_LINK_PATH TRUE
  )


# Executable to perform ray tracing
add_executable(rti_m rti/main_rt.cpp)
target_include_directories(rti_m
  PRIVATE
	  ${CMAKE_CURRENT_LIST_DIR}
    ${EMBREE_INCLUDE_DIRS}
	  ${GMSH_INCLUDE_DIR}
    ${TBB_INCLUDE_DIR}
  )
target_link_libraries(rti_m
  PRIVATE
    Boost::log
    Boost::system
    Boost::timer
    ${EMBREE_LIBRARY}
		${GMSH_LIBRARY}
    ${TBB_LIBRARY}
  )
# Copy executable straight into the directory given by the variable CMAKE_INSTALL_PREFIX.
install(TARGETS rti_m
  DESTINATION "."
  # If one wants to install the executable somewhere else, then one should consider
  # to include GNUInstallDirs and the variable CMAKE_INSTALL_BINDIR.
  )
# One needs to set the RPATH, otherwise the installed executable does not have dynamic
# link libraries set correctly.
set_target_properties(rti_m
  PROPERTIES
  INSTALL_RPATH_USE_LINK_PATH TRUE
  )







#message(STATUS "[src/CMakeLists.txt] CMAKE_INSTALL_PREFIX == ${CMAKE_INSTALL_PREFIX}")
#message(STATUS "[src/CMakeLists.txt] CMAKE_INCLUDE_PATH == ${CMAKE_INCLUDE_PATH}")
#message(STATUS "[src/CMakeLists.txt] CMAKE_LIBRARY_PATH == ${CMAKE_LIBRARY_PATH}")
#message(STATUS "[src/CMakeLists.txt] BOOST_ROOT == ${BOOST_ROOT}")
#message(STATUS "[src/CMakeLists.txt] BOOST_COMPONENTS_REQUIRED == ${BOOST_COMPONENTS_REQUIRED}")
#message(STATUS "[src/CMakeLists.txt] EMBREE_DIR == ${EMBREE_DIR}")
#message(STATUS "[src/CMakeLists.txt] EMBREE_INCLUDE_DIRS == ${EMBREE_INCLUDE_DIRS}")
#message(STATUS "[src/CMakeLists.txt] EMBREE_LIBRARY == ${EMBREE_LIBRARY}")
#message(STATUS "[src/CMakeLists.txt] EMBREE_LIBRARIES == ${EMBREE_LIBRARIES}")
#message(STATUS "[src/CMakeLists.txt] GMSH_DIR == ${GMSH_DIR}")
#message(STATUS "[src/CMakeLists.txt] GMSH_INCLUDEDIR == ${GMSH_INCLUDEDIR}")
#message(STATUS "[src/CMakeLists.txt] GMSH_LIBRARYDIR == ${GMSH_LIBRARYDIR}")
#message(STATUS "[src/CMakeLists.txt] GMSH_INCLUDE_DIR == ${GMSH_INCLUDE_DIR}")
#message(STATUS "[src/CMakeLists.txt] GMSH_LIBRARY == ${GMSH_LIBRARY}")
#message(STATUS "[src/CMakeLists.txt] TBB_DIR == ${TBB_DIR}")
#message(STATUS "[src/CMakeLists.txt] TBB_LIBRARY == ${TBB_LIBRARY}")
#message(STATUS "[src/CMakeLists.txt] TBB_LIBRARIES == ${TBB_LIBRARIES}")
#message(STATUS "[src/CMakeLists.txt] TBB_FOUND == ${TBB_FOUND}")
#message(STATUS "[src/CMakeLists.txt] CUSTOM_TBB_LIBRARY == ${CUSTOM_TBB_LIBRARY}")
