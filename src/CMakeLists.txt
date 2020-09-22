cmake_minimum_required (
  VERSION 3.11 FATAL_ERROR
  )

project (
  rti
  VERSION 0.1
  LANGUAGES CXX
  )
include(CMakePackageConfigHelpers)
include(GNUInstallDirs) # provides CMAKE_INSTALL_<dir> variables

find_package (
  embree 3 REQUIRED
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
  
# For using the library with the development executable one needs to set the RPATH,
# otherwise the installed executable does not have dynamic link libraries set correctly.
# CMAKE_INSTALL_RPATH_USE_LINK_PATH needs to be set before the executables are added.
set (
  CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE
  )

add_library (
  rtidevice INTERFACE)
target_include_directories (
  rtidevice
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  ${EMBREE_INCLUDE_DIRS}
  ${VTK_INCLUDE_DIRS}
  )
target_link_libraries (
  rtidevice
  INTERFACE
  OpenMP::OpenMP_CXX
  ${EMBREE_LIBRARIES}
  ${VTK_LIBRARIES}
  )
target_compile_features (
  rtidevice
  INTERFACE
  cxx_std_11
  )
install (
  TARGETS rtidevice
  EXPORT rtidevice-export
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
install (
  EXPORT rtidevice-export
  FILE rtideviceTargets.cmake
  NAMESPACE rtidevice::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rtidevice
  )
configure_package_config_file (
  ${CMAKE_SOURCE_DIR}/../cmake/rtideviceConfig.cmake.in
  ${CMAKE_BINARY_DIR}/cmake/rtideviceConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rtidevice
  )
install (
  FILES ${CMAKE_BINARY_DIR}/cmake/rtideviceConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rtidevice
  )
install (
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.hpp"
  )

##################
### Dev Executable
##################
add_executable (
  rti "rti/main_rt.cpp"
  )
target_link_libraries (
  rti
  PRIVATE
  rtidevice
  )
install (
  TARGETS rti
  RUNTIME
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
