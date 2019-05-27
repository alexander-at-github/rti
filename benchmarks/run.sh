#!/bin/bash

# Echo each command
set -x

SPWD=$(pwd)
MSH_FILE_NAME=cylinder.msh
GIT_HASH=$(git rev-parse HEAD)

for dir in */; do
  if [ "$dir" = "tmp/" ]; then
    continue
  fi
  cd "${SPWD}"
  cd "${dir}"
  for thrds in $(seq 1 8); do
    LOG_FILE_NAME="log.thrds.${thrds}.txt"
    date > "${LOG_FILE_NAME}"
    echo "Git HEAD hash: ${GIT_HASH}" >> "${LOG_FILE_NAME}"
    echo "Mesh file name: ${MSH_FILE_NAME}" >> "${LOG_FILE_NAME}"
    echo "Maximum number of threads: ${thrds}" >> "${LOG_FILE_NAME}"
    ./rtiex --msh-file ${MSH_FILE_NAME} --max-threads "${thrds}" >> "${LOG_FILE_NAME}"
  done
done

