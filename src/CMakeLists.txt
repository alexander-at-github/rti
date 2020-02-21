cmake_minimum_required (
  VERSION 3.11 FATAL_ERROR
  )

project (
  rti
  VERSION 0.1
  LANGUAGES CXX
  )

set (
  BUILD_SHARED_LIBS ON # A particularity of Boost 1.70
  )

find_package (
  Boost REQUIRED COMPONENTS system # "${BOOST_COMPONENTS_REQUIRED}"
  PATHS ${BOOST_ROOT}
  NO_DEFAULT_PATH
  )
find_package (
  embree 3.6 REQUIRED
  PATHS ${EMBREE_DIR}
  NO_DEFAULT_PATH
  )
find_package (
  VTK 8.2 REQUIRED
  PATHS ${VTK_DIR}
  NO_DEFAULT_PATH
  )
find_package (
  OpenMP REQUIRED
  )

# One needs to set the RPATH, otherwise the installed executable does not have dynamic
# link libraries set correctly.
# CMAKE_INSTALL_RPATH_USE_LINK_PATH needs to be set before the executables are added.
set (
  CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE
  )

### Create library
add_library (
  librtidevice STATIC "rti/device.cpp"
  )
set_target_properties (
  librtidevice
  PROPERTIES PUBLIC_HEADER "rti/device.hpp"
  )
target_include_directories (
  librtidevice
  PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}
  PUBLIC
  ${EMBREE_INCLUDE_DIRS}
  ${VTK_INCLUDE_DIRS}
  $<INSTALL_INTERFACE:include>
  )
target_compile_features (
  librtidevice
  PRIVATE
  cxx_std_11
  )
target_link_libraries (
  librtidevice
  PUBLIC
  OpenMP::OpenMP_CXX
  Boost::system
  ${EMBREE_LIBRARIES}
  ${VTK_LIBRARIES}
  )

### Export the library
include(GNUInstallDirs)
# GNUInstallDirs provides CMAKE_INSTALL_<dir> variables
install (
  TARGETS librtidevice
  EXPORT librtidevice-export
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rti
  )
install (
  EXPORT librtidevice-export
  FILE librtideviceTargets.cmake
  NAMESPACE librtidevice::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/librtidevice
  )

### Create a Config file for the library for other CMake projects to use find_package()
include(CMakePackageConfigHelpers)
configure_package_config_file (
  ${CMAKE_SOURCE_DIR}/../cmake/librtideviceConfig.cmake.in
  ${CMAKE_BINARY_DIR}/cmake/librtideviceConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/librtidevice
  )
install (
  FILES ${CMAKE_BINARY_DIR}/cmake/librtideviceConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/librtidevice
  )
