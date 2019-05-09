externalproject_add(tbb_external
  # Only download
  GIT_REPOSITORY
    https://github.com/01org/tbb.git
  # Hash of tags/2019_U5
  GIT_TAG
    112ade3da541b73b053635b6f29cafe4abe4fa5e
  # No configure
  CONFIGURE_COMMAND ""
  # No build
  BUILD_COMMAND ""
  # No install
  INSTALL_COMMAND
    ""
  # TBB does not need to be build. We use it only with embree and embree builds TBB itself!
  #CONFIGURE_COMMAND
  #  ${CMAKE_COMMAND} -E "copy" ${CMAKE_CURRENT_LIST_DIR}/TBBCMakeLists.txt <SOURCE_DIR>/CMakeLists.txt
  #  COMMAND ${CMAKE_COMMAND} -B . -S <SOURCE_DIR>

  #CMAKE_ARGS
  #  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  #  -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
  #  -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
  #  -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
  #  #-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> # strange that we need that at all
  #CMAKE_CACHE_ARGS
  #  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  #  -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
  USES_TERMINAL_DOWNLOAD
    1
  USES_TERMINAL_UPDATE
    1
  GIT_PROGRESS
    1
)

#externalproject_get_property(tbb_external SOURCE_DIR BINARY_DIR)
#set(TBB_SOURCE_DIR ${SOURCE_DIR})
#set(TBB_BINARY_DIR ${BINARY_DIR})
#unset(SOURCE_DIR)
#unset(BINARY_DIR)
#
#message(STATUS "[tbb/CMakeLists.txt] TBB_BINARY_DIR == ${TBB_BINARY_DIR}")
#message(STATUS "[tbb/CMakeLists.txt] TBB_SOURCE_DIR == ${TBB_SOURCE_DIR}")
#
#set(
#  TBB_DIR "${TBB_SOURCE_DIR}/cmake"
#  CACHE PATH "Path to TBBConfig.cmake"
#  FORCE
#  )
#set(
#  TBB_LIBRARYDIR "${TBB_BINARY_DIR}/tbb_cmake_build/tbb_cmake_build_subdir_release"
#  CACHE PATH "Path to TBB library directories"
#  FORCE
#  )
#
#message(STATUS "[tbb/CMakeLists.txt] TBB_LIBRARYDIR == ${TBB_LIBRARYDIR}")
#message(STATUS "[tbb/CMakeLists.txt] TBB_DIR == ${TBB_DIR}")
