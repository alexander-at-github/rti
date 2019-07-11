# Get OCE installation directory for Gmsh
externalproject_get_property(oce-external INSTALL_DIR)
# rename the variable to something more meaningful
set(OCE_INSTALL_PREFIX ${INSTALL_DIR})
unset(INSTALL_DIR)

externalproject_add(gmsh-external
  DEPENDS
    oce-external
  GIT_REPOSITORY
    https://gitlab.onelab.info/gmsh/gmsh.git
  # Hash of tags/gmsh_4_3_0
  GIT_TAG
    8375000a0c4b3d9b54e7d1c39fe6b104b9c34a5d
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    -DENABLE_FLTK=1 # One may disable the compilation of the GUI components, if not needed.
    #-DENABLE_BUILD_LIB=1
    -DENABLE_BUILD_SHARED=1
    # Enable dynamic Gmsh executable (linked with shared lib) (default: OFF)
    -DENABLE_BUILD_DYNAMIC=1
  CMAKE_CACHE_ARGS
    #-DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
    -DENV{CASROOT}:PATH=${OCE_INSTALL_PREFIX}
  USES_TERMINAL_DOWNLOAD
    1
  USES_TERMINAL_UPDATE
    1
  GIT_PROGRESS
    1
)

# Get install and build directories from Gmsh
externalproject_get_property(gmsh-external INSTALL_DIR BINARY_DIR)

set(
  GMSH_DIR "${INSTALL_DIR}"
  CACHE PATH "Path to Gmsh installation"
  FORCE
  )
#set(
#  GMSH_INCLUDEDIR "${GMSH_DIR}/include"
#  CACHE PATH "Path to Gmsh include directories"
#  FORCE
#  )
#set(
#  GMSH_LIBRARYDIR "${GMSH_DIR}/lib"
#  CACHE PATH "Path to Gmsh library directories"
#  FORCE
#  )
