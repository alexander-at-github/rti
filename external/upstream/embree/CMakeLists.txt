# Embree-TBB specific
externalproject_get_property(tbb-external SOURCE_DIR)
set(TBB_SOURCE_DIR ${SOURCE_DIR})
unset(SOURCE_DIR)

#message(STATUS "[embree/CMakeLists.txt] TBB_SOURCE_DIR=${TBB_SOURCE_DIR}")

externalproject_add(embree-external
  DEPENDS
    tbb-external
  GIT_REPOSITORY
    https://github.com/embree/embree.git
  # Hash of tags/v3.5.2
  GIT_TAG
    e766297adcdaf8e5d8dd996c15aed753255230fd
  CMAKE_ARGS
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> # strange that we need that at all
    -DEMBREE_TUTORIALS=OFF

    #-DCMAKE_BUILD_TYPE=Debug
  CMAKE_CACHE_ARGS
    #-DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
    -DEMBREE_ISPC_SUPPORT:BOOL=${EMBREE_ISPC_SUPPORT}
    # Embree TBB specific
    # Embree reads this variable in its own CMakeLists.txt.
    # (TBB_ROOT is not used according to the convetion of CMake).
    -DTBB_ROOT:PATH=${TBB_SOURCE_DIR}
  # Do not install
  # If one installes it and tries to use the installed library, then there exists a linkage
  # problem with TBB. Simply use the generated files in the binary directory.
  INSTALL_COMMAND
    ""
  USES_TERMINAL_DOWNLOAD
    1
  USES_TERMINAL_UPDATE
    1
  GIT_PROGRESS
    1
)

# Get install and build directory from external project
externalproject_get_property(embree-external INSTALL_DIR BINARY_DIR)
#message(STATUS "[embree/CMakeLists.txt] INSTALL_DIR=${INSTALL_DIR}")
#message(STATUS "[embree/CMakeLists.txt] BINARY_DIR=${BINARY_DIR}")

set(
  # Very important to use the binary directory as EMBREE_DIR here, cause when the Embree
  # library is copied to the install directory its dynamic link to the TBB library seems
  # to break.
  EMBREE_DIR "${BINARY_DIR}"
  CACHE PATH "Path to Embree installation"
  FORCE
  )
set(
  EMBREE_INCLUDEDIR "${EMBREE_DIR}/include"
  CACHE PATH "Path to Embree include directories"
  FORCE
  )
set(
  EMBREE_LIBRARYDI "${EMBREE_DIR}/lib"
  CACHE PATH "Path to Embree library directories"
  FORCE
  )
