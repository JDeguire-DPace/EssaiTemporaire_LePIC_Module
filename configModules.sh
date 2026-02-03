#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="${1:-build}"
BUILD_TYPE="${2:-Release}"

# If your cluster requires modules, uncomment and adapt:
# module load intel-oneapi
# module load intel-oneapi-mpi

cmake -S . -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DCMAKE_Fortran_COMPILER=mpiifx \
  -DENABLE_OPENMP=ON

cmake --build "${BUILD_DIR}" -j
echo "Built: ${BUILD_DIR}/LePIC_3D"