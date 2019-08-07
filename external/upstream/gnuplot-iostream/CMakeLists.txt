externalproject_add(gnuplot-iostream-external
  GIT_REPOSITORY
  https://github.com/dstahlke/gnuplot-iostream.git
  GIT_TAG
  # Latest commit at the time of writing
  7a5fed24bbce977caafeeed773b20dbf31f43d52
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

externalproject_get_property(gnuplot-iostream-external SOURCE_DIR)
set(
  GPLOT_IO_DIR "${SOURCE_DIR}"
  CACHE PATH "Path to gnuplot-iostream source"
  FORCE
  )
