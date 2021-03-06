cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(
  rti-bench
  VERSION 1.0
  LANGUAGES C CXX
  )

#message(STATUS "[test/CMakeLists.txt] GBENCH_SRC_DIR == ${GBENCH_SRC_DIR}")
#message(STATUS "[test/CMakeLists.txt] GBENCH_BIN_DIR == ${GBENCH_BIN_DIR}")
#message(STATUS "[test/CMakeLists.txt] GBENCH_INST_DIR == ${GBENCH_INST_DIR}")
#message(STATUS "[test/CMakeLists.txt] BOOST_ROOT == ${BOOST_ROOT}")

# One needs to set the RPATH, otherwise the installed executable does not have dynamic
# link libraries set correctly.
# CMAKE_INSTALL_RPATH_USE_LINK_PATH needs to be set before the executables are added.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(benchmark REQUIRED
  PATHS ${GBENCH_INST_DIR}
  NO_DEFAULT_PATH
  )
find_package(embree 3 REQUIRED
  PATHS ${EMBREE_DIR}
  NO_DEFAULT_PATH
  )


add_executable(benchmark "")
target_sources(benchmark
  PRIVATE
  rti/dummy_benchmark.cpp
  rti/ray/rectangle_origin_z.cpp
  rti/intersect_vs_occluded_all.cpp
  )
target_include_directories(benchmark
  PRIVATE
  ${RTI_SRC_DIR}
  )
target_link_libraries(benchmark
  PRIVATE
  benchmark::benchmark
  # also takes care of include directories
  benchmark::benchmark_main
  #
  ${EMBREE_LIBRARIES}
  )
install(
  TARGETS benchmark
  DESTINATION ".")
