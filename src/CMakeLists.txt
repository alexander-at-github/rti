cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(
  rti
  VERSION 0.1
  LANGUAGES CXX
  )

# message(STATUS "CMAKE_CXX_FLAGS == ${CMAKE_CXX_FLAGS}")
# message(STATUS "CMAKE_BUILD_TYPE == '${CMAKE_BUILD_TYPE}'")
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#   message(STATUS "CMAKE_CXX_FLAGS_DEBUG == ${CMAKE_CXX_FLAGS_DEBUG}")
# endif()
# if(CMAKE_BUILD_TYPE STREQUAL "Release")
#   message(STATUS "CMAKE_CXX_FLAGS_RELEASE == ${CMAKE_CXX_FLAGS_RELEASE}")
# endif()

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
find_path(GMSH_INCLUDE_DIRS NAMES gmsh.h
  PATHS "${GMSH_DIR}/include"
	NO_DEFAULT_PATH
	)
find_library(GMSH_LIBRARIES NAMES gmsh libgmsh
  PATHS "${GMSH_DIR}/lib"
  NO_DEFAULT_PATH
  )
# message(STATUS "VTK_DIR == ${VTK_DIR}")
find_package(VTK 8.2 REQUIRED
  PATHS ${VTK_DIR}
  NO_DEFAULT_PATH
  )

find_package(OpenMP REQUIRED)

# One needs to set the RPATH, otherwise the installed executable does not have dynamic
# link libraries set correctly.
# CMAKE_INSTALL_RPATH_USE_LINK_PATH needs to be set before the executables are added.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# # Executable for testing the ray tracing code
# add_executable(rtiex rti/main_testing_rt.cpp)
# target_include_directories(rtiex
#   PRIVATE
# 	  ${CMAKE_CURRENT_LIST_DIR}
#     ${EMBREE_INCLUDE_DIRS}
# 	  ${GMSH_INCLUDE_DIR}
#   )
# target_link_libraries(rtiex
#   PRIVATE
#     OpenMP::OpenMP_CXX # this line also takes care for flags
#     Boost::log
#     Boost::system
#     Boost::timer
#     ${EMBREE_LIBRARY}
# 		${GMSH_LIBRARY}
#   )
# # Copy executable straight into the directory given by the variable CMAKE_INSTALL_PREFIX.
# install(TARGETS rtiex
#   DESTINATION "."
#   # If one wants to install the executable somewhere else, then one should consider
#   # to include GNUInstallDirs and the variable CMAKE_INSTALL_BINDIR.
#   )
# # One needs to set the RPATH, otherwise the installed executable does not have dynamic
# # link libraries set correctly.
# set_target_properties(rtiex
#   PROPERTIES
#   INSTALL_RPATH_USE_LINK_PATH TRUE
#   )


# Executable to perform ray tracing
add_executable(rti_m rti/main_rt.cpp)
target_include_directories(rti_m
  PRIVATE
	  ${CMAKE_CURRENT_LIST_DIR}
    ${EMBREE_INCLUDE_DIRS}
	  ${GMSH_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
  )
target_link_libraries(rti_m
  PRIVATE
    OpenMP::OpenMP_CXX # this line also takes care for flags
    Boost::log
    Boost::system
    Boost::timer
    ${EMBREE_LIBRARIES}
		${GMSH_LIBRARIES}
    ${VTK_LIBRARIES}
  )
# # One needs to set the RPATH, otherwise the installed executable does not have dynamic
# # link libraries set correctly.
# set_target_properties(rti_m
#   PROPERTIES
#   INSTALL_RPATH_USE_LINK_PATH TRUE
#   )



# Executable to read text input files
include(${VTK_USE_FILE}) # VTK specific
add_executable(read_dsv rti/dsv_reader.cpp)
target_include_directories(read_dsv
  PRIVATE
  ${VTK_INCLUDE_DIRS}
  )
target_link_libraries(read_dsv
  PRIVATE
  ${VTK_LIBRARIES}
  )



# Copy executable straight into the directory given by the variable CMAKE_INSTALL_PREFIX.
install(TARGETS rti_m read_dsv
  DESTINATION "."
  # If one wants to install the executable somewhere else, then one should consider
  # to include GNUInstallDirs and the variable CMAKE_INSTALL_BINDIR.
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
