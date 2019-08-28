cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(
  rti-benchmark
  VERSION 1.0
  LANGUAGES CXX
  )

# One needs to set the RPATH, otherwise the installed executable does not have dynamic
# link libraries set correctly.
# CMAKE_INSTALL_RPATH_USE_LINK_PATH needs to be set before the executables are added.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(embree 3 REQUIRED
  PATHS ${EMBREE_DIR}
  NO_DEFAULT_PATH
  )


add_executable(multi-intersect rti/multi_intersect.cpp)
target_include_directories(multi-intersect
	PRIVATE
  ${RTI_SRC_DIR}
  )
target_link_libraries(multi-intersect
  PRIVATE
  ${EMBREE_LIBRARIES}
  )
add_executable(intersect-all rti/intersect_all.cpp)
target_include_directories(intersect-all
  PRIVATE
  ${RTI_SRC_DIR}
  )
target_link_libraries(intersect-all
  PRIVATE
  ${EMBREE_LIBRARIES}
  )

# Copy executable straight into the directory given by the variable CMAKE_INSTALL_PREFIX.
install(TARGETS multi-intersect intersect-all
  DESTINATION "."
  # If one wants to install the executable somewhere else, then one should consider
  # to include GNUInstallDirs and the variable CMAKE_INSTALL_BINDIR.
  )

