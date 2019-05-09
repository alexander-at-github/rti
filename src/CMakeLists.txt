cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(
  rti_core
  VERSION 0.1
  LANGUAGES CXX
  )

set(BUILD_SHARED_LIBS ON) # A particularity of Boost 1.70

find_package(Boost REQUIRED COMPONENTS log system # "${BOOST_COMPONENTS_REQUIRED}"
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

add_executable(rtiex rti/main.cpp)

target_include_directories(rtiex
  PRIVATE
	  ${CMAKE_CURRENT_LIST_DIR}
	  ${GMSH_INCLUDE_DIR}
    ${EMBREE_INCLUDE_DIRS}
  )
target_link_libraries(rtiex
  PRIVATE
    Boost::log
    Boost::system
		${GMSH_LIBRARY}
    ${EMBREE_LIBRARY}
  )

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
