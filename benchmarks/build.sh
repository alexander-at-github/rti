#!/bin/bash

# Echo each command
set -x

SPWD=$(pwd)

for dir in */; do
  if [ "$dir" = "tmp/" ]; then
    continue
  fi
  cd "${SPWD}"
  cd "${dir}"
  rm -r build.ninja CMakeCache.txt CMakeFiles/ cmake_install.cmake external/ rtiex rti-prefix/ rules.ninja .ninja*
  cmake -DCMAKE_BUILD_TYPE=Release -G Ninja ../..
  cmake --build . --target rti
done
