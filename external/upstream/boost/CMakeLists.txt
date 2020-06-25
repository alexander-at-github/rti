### Boost does not use Cmake
### We need to use a custom build

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  message(STATUS "boost cmake build recognized GNU compiler")
  if(APPLE)
    set(_toolset "darwin")
  else()
    set(_toolset "gcc")
  endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  message(STATUS "boost cmake build recognized Clang compiler")
  set(_toolset "clang")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  message(STATUS "boost cmake build recognized Intel compiler")
  if(APPLE)
    set(_toolset "intel-darwin")
  else()
    set(_toolset "intel-linux")
  endif()
endif()

if(NOT "${BOOST_COMPONENTS_REQUIRED}" STREQUAL "")
  # Replace unit_test_framework (used by CMake's find_package) with test (understood by Boost build toolchain)
  string(REPLACE "unit_test_framework" "test" _b2_needed_components "${BOOST_COMPONENTS_REQUIRED}")
  # Generate argument for BUILD_BYPRODUCTS
  set(_build_byproducts)
  set(_b2_select_libraries)
  foreach(_lib IN LISTS _b2_needed_components)
    list(APPEND _build_byproducts ${STAGED_INSTALL_PREFIX}/boost/lib/libboost_${_lib}${CMAKE_SHARED_LIBRARY_SUFFIX})
    list(APPEND _b2_select_libraries --with-${_lib})
  endforeach()
  # Transform the ;-separated list to a ,-separated list (digested by the Boost build toolchain!)
  string(REPLACE ";" "," _b2_needed_components "${_b2_needed_components}")
  set(_bootstrap_select_libraries "--with-libraries=${_b2_needed_components}")
  string(REPLACE ";" ", " printout "${BOOST_COMPONENTS_REQUIRED}")
  message(STATUS "  Libraries to be built: ${printout}")
endif()

ExternalProject_Add(boost-external
  URL
    #https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.zip
    https://sourceforge.net/projects/boost/files/boost/1.70.0/boost_1_70_0.tar.gz
  URL_HASH
    #SHA256=02d420e6908016d4ac74dfc712eec7d9616a7fc0da78b0a1b5b937536b2e01e8
    SHA256=882b48708d211a5f48e60b0124cf5863c1534cd544ecd0664bb534a4b5d506e9
  #DOWNLOAD_NO_PROGRESS
  #  1
  UPDATE_COMMAND
    ""
  CONFIGURE_COMMAND
    <SOURCE_DIR>/bootstrap.sh
    --with-toolset=${_toolset}
    --prefix=${STAGED_INSTALL_PREFIX}/boost
    ${_bootstrap_select_libraries}
  BUILD_COMMAND
    <SOURCE_DIR>/b2 -q
    link=shared
    threading=multi
    variant=release
    toolset=${_toolset}
    ${_b2_select_libraries}
    -j 8
  LOG_BUILD
    1
  BUILD_IN_SOURCE
    1
  INSTALL_COMMAND
    <SOURCE_DIR>/b2 -q install
    link=shared
    threading=multi
    variant=release
    toolset=${_toolset}
    ${_b2_select_libraries}
  LOG_INSTALL
    1
  BUILD_BYPRODUCTS
    "${_build_byproducts}"
  )

set(
  BOOST_ROOT ${STAGED_INSTALL_PREFIX}/boost
  CACHE PATH "Path to internally built Boost installation root"
  FORCE
  )
set(
  BOOST_INCLUDEDIR ${BOOST_ROOT}/include
  CACHE PATH "Path to internally built Boost include directories"
  FORCE
  )
set(
  BOOST_LIBRARYDIR ${BOOST_ROOT}/lib
  CACHE PATH "Path to internally built Boost library directories"
  FORCE
  )

unset(_toolset)
unset(_b2_needed_components)
unset(_build_byproducts)
unset(_b2_select_libraries)
unset(_boostrap_select_libraries)
