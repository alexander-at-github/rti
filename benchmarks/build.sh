#!/bin/bash

# Echo each command
set -x

SPWD=$(pwd)

for dir in ./do/*/; do
  #if [ "$dir" = "tmp/" ]; then
  #  continue
  #fi
  cd "${SPWD}"
  cd "${dir}"
  rm -r build.ninja CMakeCache.txt CMakeFiles/ cmake_install.cmake external/ rtiex rti-prefix/ rules.ninja .ninja*
  cmake -DCMAKE_BUILD_TYPE=Release -G Ninja ../../..
  #cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-4.9 -DCMAKE_CXX_COMPILER=g++-4.9 -G Ninja ../../..
  #cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=/home/alexander/local/opt/intel/system_studio_2019/bin/icc -DCMAKE_CXX_COMPILER=/home/alexander/local/opt/intel/system_studio_2019/bin/icc -G Ninja ../../..
  cmake --build . --target rti
done
