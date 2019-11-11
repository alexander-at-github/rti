externalproject_add(embree-external
  GIT_REPOSITORY
    https://github.com/embree/embree.git
  GIT_TAG
    ## Hash of tags/v3.5.2
    #e766297adcdaf8e5d8dd996c15aed753255230fd
    #
    # Hash of tags/v3.6.1
    d17aec1491a6c664f4da2f0b8ba03c171cbf36a1
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    -DCMAKE_SHARED_LINKER_FLAGS=-lstdc++
    -DEMBREE_TUTORIALS=OFF
    -DEMBREE_BACKFACE_CULLING=ON
    -DEMBREE_FILTER_FUNCTION=ON # setting this option is probably not necessary

    #-DCMAKE_BUILD_TYPE=Debug
  CMAKE_CACHE_ARGS
    #-DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
    -DEMBREE_ISPC_SUPPORT:BOOL=${EMBREE_ISPC_SUPPORT}
    # We use OpenMP, Embree calls that internal (eventhough it is external to Embree)
    -DEMBREE_TASKING_SYSTEM:STRING=INTERNAL
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

# The following might be simplified, because we are not using TBB anymore.

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
