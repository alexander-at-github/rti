cmake_minimum_required(VERSION 3.7)
project(rti)

## Add libraries
add_subdirectory(lib)
## Add package Embree 3.5
set(embree_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/embree")
find_package(embree 3.5.2 REQUIRED)
## Add Embree headers to include directories
include_directories("${EMBREE_INCLUDE_DIRS}")

file(GLOB SOURCES "src/**/*.cpp")
add_executable(TestRTIntersect ${SOURCES})
