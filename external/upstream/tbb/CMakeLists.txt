## Convert compiler ID to TBB format.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #MESSAGE("Clang")
  set(TBB_COMPILER "clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  #MESSAGE("Clang")
  set(TBB_COMPILER "gcc")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  #MESSAGE("Clang")
  set(TBB_COMPILER "icc")
endif()

# Set TBB build options
set(TBB_COMPONENTS tbb tbbmalloc)
#set(TBB_BUILD_TYPE debug)
#set(TBB_BUILD_TYPE release)

externalproject_add(tbb-external
  # Only download
  GIT_REPOSITORY
    https://github.com/01org/tbb.git
  # Hash of tags/2019_U5
  GIT_TAG
  112ade3da541b73b053635b6f29cafe4abe4fa5e
  ## TBB Build directions are given in <tbb-source-dir>/build/index.html
  # No configure
  CONFIGURE_COMMAND
    ""
  #
  BUILD_COMMAND
  cd <SOURCE_DIR>
  #COMMAND make ${TBB_COMPONENTS} compiler=${TBB_COMPILER} stdver=c++${CMAKE_CXX_STANDARD} tbb_build_dir=<BINARY_DIR>
  COMMAND make ${TBB_COMPONENTS} compiler=${TBB_COMPILER} stdver=c++${CMAKE_CXX_STANDARD} tbb_build_prefix=lib

  # No install
  INSTALL_COMMAND
    ""

  # Note: The build system of Embree also builds TBB for its own use.

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

externalproject_get_property(tbb-external SOURCE_DIR)
set(
  TBB_DIR ${SOURCE_DIR}
  CACHE PATH "TBB (source) directory"
  FORCE
  )
unset(SOURCE_DIR)
