cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

set (NUM_TRUNC_CHARS 2)
set (RPATH "R")
set (RSCRIPT_PATH "Rscript")
set (RREPO "https://cran.wu.ac.at/")
# execute_process (COMMAND ${RSCRIPT_PATH} -e "install.packages('Rcpp', repos='${RREPO}')")
# execute_process (COMMAND ${RSCRIPT_PATH} -e "install.packages('RcppArmadillo', repos='${RREPO}')")
# execute_process (COMMAND ${RSCRIPT_PATH} -e "install.packages('RInside', repos='${RREPOS}')")

set (INCLUDE_DIRS "") # empty list


execute_process(
  COMMAND ${RPATH} CMD config --cppflags
  OUTPUT_VARIABLE RINCL)
string(STRIP ${RINCL} RINCL)
string(SUBSTRING ${RINCL} ${NUM_TRUNC_CHARS} -1 RINCL)
list(APPEND INCLUDE_DIRS ${RINCL})
message(STATUS "RINCL == ${RINCL}")

execute_process(
  COMMAND ${RSCRIPT_PATH} -e "Rcpp:::CxxFlags()"
  OUTPUT_VARIABLE RCPPINCL)
string(SUBSTRING ${RCPPINCL} ${NUM_TRUNC_CHARS} -1 RCPPINCL)
list(APPEND INCLUDE_DIRS ${RCPPINCL})

execute_process(
  COMMAND ${RSCRIPT_PATH} -e "RInside:::CxxFlags()"
  OUTPUT_VARIABLE RINSIDEINCL)
string(SUBSTRING ${RINSIDEINCL} ${NUM_TRUNC_CHARS} -1 RINSIDEINCL)
list(APPEND INCLUDE_DIRS ${RINSIDEINCL})

execute_process(
  COMMAND ${RSCRIPT_PATH} -e "RcppArmadillo:::CxxFlags()"
  OUTPUT_VARIABLE RCPPARMADILLOINCL)
message(STATUS "RCPPARMADILLOINCL ==${RCPPARMADILLOINCL}")
string(LENGTH ${RCPPARMADILLOINCL} INCLLENGTH)
math(EXPR INCLLENGTH "${INCLLENGTH}-4")
string(SUBSTRING ${RCPPARMADILLOINCL} 3 ${INCLLENGTH} RCPPARMADILLOINCL)
list(APPEND INCLUDE_DIRS ${RCPPARMADILLOINCL})

set (LIBRARIES "") # empty list

execute_process(
  COMMAND ${RPATH} CMD config --ldflags
  OUTPUT_VARIABLE RLIBS)
string(STRIP ${RLIBS} RLIBS)
message(STATUS "RLIBS == ${RLIBS}")
string(LENGTH ${RLIBS} RLIBS_LEN)
if (${RLIBS} MATCHES "[-][L]([^ ;])+")
  string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLIBS_L)
  string(STRIP ${RLIBS_L} RLIBS_L )
endif()
if (${RLIBS} MATCHES "[-][l]([^;])+")
  string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLIBS_l)
  string(STRIP ${RLIBS_l} RLIBS_l )
endif()
list(APPEND LIBRARIES ${RLIBS})

execute_process(
  COMMAND ${RSCRIPT_PATH} -e "Rcpp:::LdFlags()"
  OUTPUT_VARIABLE RCPPLIBS)
if (${RCPPLIBS} MATCHES "[-][L]([^ ;])+")
  string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPPLIBS_L)
  link_directories(${RCPPLIBS_L} )
endif()
if (${RCPPLIBS} MATCHES "[-][l][R]([^;])+")
  string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPPLIBS_l)
endif()
message(STATUS "RCPPLIBS_L == ${RCPPLIBS_L}")
message(STATUS "RCPPLIBS_l == ${RCPPLIBS_l}")
list(APPEND LIBRARIES ${RCPPLIBS})
execute_process(
  COMMAND ${RSCRIPT_PATH} -e "RInside:::LdFlags()"
  OUTPUT_VARIABLE RINSIDELIBS)
if (${RINSIDELIBS} MATCHES "[-][L]([^ ;])+")
  string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDELIBS_L)
  link_directories(${RINSIDELIBS_L})
endif()
if (${RINSIDELIBS} MATCHES "[-][l][R]([^;])+")
  string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDELIBS_l)
endif()
message(STATUS "RINSIDELIBS_L == ${RINSIDELIBS_L}")
message(STATUS "RINSIDELIBS_l == ${RINSIDELIBS_l}")
list(APPEND LIBRARIES ${RINSIDELIBS})

set (
  RINSIDE_INCLUDE_DIRS ${INCLUDE_DIRS}
  CACHE PATH "Paths to RInside include directories"
  FORCE
  )
message(STATUS "RINSIDE_INCLUDE_DIRS == ${RINSIDE_INCLUDE_DIRS}")

set (
  RINSIDE_LIBRARIES ${LIBRARIES}
  CACHE PATH "Paths to RInside library locations"
  FORCE
  )
message(STATUS "RINSIDE_LIBRARIES == ${RINSIDE_LIBRARIES}")
