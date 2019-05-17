## For building gperftools one needs autoconf, automake, and libtool

externalproject_add(gperftools-external
  GIT_REPOSITORY
  https://github.com/gperftools/gperftools.git
  # Hash of tags/gperftools-2.7
  GIT_TAG
  9608fa3bcf8020d35f59fbf70cd3cbe4b015b972
  #
  CONFIGURE_COMMAND
  cd <SOURCE_DIR>
  COMMAND autoreconf -fvi # hack to fix errors
  COMMAND ./configure --prefix=<INSTALL_DIR>
  #
  BUILD_COMMAND
  cd <SOURCE_DIR>
  COMMAND make
  #
  INSTALL_COMMAND
  cd <SOURCE_DIR>
  COMMAND make install
  # For more information on building gperftools see the file
  # INSTALL in the gperftools source directory.
  )
