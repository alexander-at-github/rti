externalproject_add(gtest-external
  GIT_REPOSITORY
  https://github.com/google/googletest.git
  GIT_TAG
  # Google Test follows the Live at Head philosophy
  master
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
    -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
    #-DCMAKE_SHARED_LINKER_FLAGS=-lstdc++
  CMAKE_CACHE_ARGS
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
  USES_TERMINAL_DOWNLOAD
    1
  USES_TERMINAL_UPDATE
    1
  GIT_PROGRESS
    1
  )


# Get install and build directories from GTest
externalproject_get_property(gtest-external SOURCE_DIR BINARY_DIR)

set(
  GTEST_SRC_DIR "${SOURCE_DIR}"
  CACHE PATH "Path to GTest source"
  FORCE
  )

set(
  GTEST_BIN_DIR "${BINARY_DIR}"
  CACHE PATH "Path to GTest build dir"
  FORCE
  )
