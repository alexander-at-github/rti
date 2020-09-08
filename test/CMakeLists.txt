cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

project(
  rti-test
  VERSION 1.0
  LANGUAGES CXX
  )

message(STATUS "[test/CMakeLists.txt] GTEST_DIR == ${GTEST_DIR}")
message(STATUS "[test/CMakeLists.txt] GTEST_BIN_DIR == ${GTEST_BIN_DIR}")
message(STATUS "[test/CMakeLists.txt] BOOST_ROOT == ${BOOST_ROOT}")
message(STATUS "[test/CMakeLists.txt] GPLOT_IO_DIR == ${GPLOT_IO_DIR}")

# One needs to set the RPATH, otherwise the installed executable does not have dynamic
# link libraries set correctly.
# CMAKE_INSTALL_RPATH_USE_LINK_PATH needs to be set before the executables are added.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# Incorporate GTest directly into the build
# See https://stackoverflow.com/questions/40267752/cmake-ctest-test-failure-when-initializing-googletest#answer-47682495
# and https://cmake.org/cmake/help/latest/module/GoogleTest.html
# and the CMake Cookbook
add_subdirectory(${GTEST_SRC_DIR} ${GTEST_BIN_DIR})
include(GoogleTest)

# gnuplot-iostream uses boost
# set shared libs influences all libraries added after this point.
# That is, it might interfer with other libraries in unexpected ways.
set(BUILD_SHARED_LIBS ON) # for Boost 1.70
find_package(Boost REQUIRED COMPONENTS system filesystem iostreams
  PATHS ${BOOST_ROOT}
  NO_DEFAULT_PATH
  )

add_executable(tests "")
target_sources(tests
  PRIVATE
  rti/ray/rectangle_origin_z.cpp
  rti/ray/cosine_direction.cpp
  rti/geo/disc_bounding_box_intersector.cpp
  )
target_include_directories(tests
  PRIVATE
  ${RTI_SRC_DIR}
  ${GPLOT_IO_DIR})
target_link_libraries(tests
  PRIVATE
  gtest_main
  Boost::iostreams
  Boost::system
  Boost::filesystem)
enable_testing()
add_test(
  NAME gtest
  COMMAND $<TARGET_FILE:tests>
  )

install(
  TARGETS tests
  DESTINATION ".")

#add_executable(test_rectangle_origin_z rti/ray/rectangle_origin_z.cpp)
#target_include_directories(test_rectangle_origin_z)
#enable_testing()
#gtest_discover_tests(test_rectangle_origin_z)
