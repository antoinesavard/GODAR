#!/bin/bash
set -euo pipefail

# ============================================================
# Build script for GODAR and its dependencies (coretran)
# on LUMI using amdflang via the mpifort wrapper.
#
# Run from any directory that contains src/ :
#   bash src/lumi/build.sh
# ============================================================

# Locate the repository root (parent of src/)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOPDIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [ ! -d "${TOPDIR}/src" ]; then
    echo "ERROR: cannot find src/ under ${TOPDIR}" >&2
    exit 1
fi

# ------------------------------------------------------------
# 0. Ensure sources use `use mpi` (not mpi_f08)
#    build-f08-bug.sh patches to mpi_f08; this reverts it.
# ------------------------------------------------------------
for _f in "${TOPDIR}"/src/model/src/{godar,stepper,diagnostics_mpi,mpi_coms}.f90; do
    sed -i 's/use mpi_f08$/use mpi/' "$_f"
done

CORETRAN_SRC="${TOPDIR}/coretran"
CORETRAN_INSTALL="${TOPDIR}/coretran-install"
GODAR_SRC="${TOPDIR}/src"
GODAR_BUILD="${GODAR_SRC}/build"
VENV_DIR="${TOPDIR}/venv"

# ------------------------------------------------------------
# 1. Load modules
# ------------------------------------------------------------
echo "=== Loading modules ==="
module use /appl/local/containers/test-modules
ml rocm-afar-cray-mpich/24.3.0
ml netcdf_fortran
ml cray-python

echo "Fortran compiler (mpifort wraps amdflang):"
mpifort --version | head -2
echo ""
echo "Python from cray-python:"
python3 --version
echo ""

# ------------------------------------------------------------
# 2. Create Python virtual environment
# ------------------------------------------------------------
if [ ! -d "${VENV_DIR}" ]; then
    echo "=== Creating Python virtual environment ==="
    python3 -m venv "${VENV_DIR}"
fi

echo "=== Activating Python virtual environment ==="
source "${VENV_DIR}/bin/activate"
echo "Python venv active: $(which python3) — $(python3 --version)"
echo ""

# ------------------------------------------------------------
# 3. Clean previous build artifacts
# ------------------------------------------------------------
echo "=== Cleaning previous builds ==="
rm -rf "${CORETRAN_SRC}/build" "${CORETRAN_INSTALL}"
rm -rf "${GODAR_BUILD}" "${GODAR_SRC}/bin" "${GODAR_SRC}/lib" "${GODAR_SRC}/include"
echo ""

# ------------------------------------------------------------
# 4. Clone and build coretran
# ------------------------------------------------------------
if [ ! -d "${CORETRAN_SRC}" ]; then
    echo "=== Cloning coretran ==="
    git clone https://github.com/leonfoks/coretran.git "${CORETRAN_SRC}"
fi

echo "=== Building coretran ==="
mkdir -p "${CORETRAN_SRC}/build"
cd "${CORETRAN_SRC}/build"

# Patch: add binary dir for out-of-tree test subdirectory
sed -i 's|add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../test)|add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../test ${CMAKE_CURRENT_BINARY_DIR}/test)|' \
    "${CORETRAN_SRC}/src/CMakeLists.txt"

cmake \
    -DCMAKE_Fortran_COMPILER=mpifort \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DCMAKE_INSTALL_PREFIX="${CORETRAN_INSTALL}" \
    -DCMAKE_Fortran_FLAGS="-fPIC" \
    -DBUILD_SHARED_LIBS=OFF \
    ../src

make -j"$(nproc)" install

echo ""
echo "=== coretran installed to ${CORETRAN_INSTALL} ==="
echo ""

# ------------------------------------------------------------
# 5. Build GODAR
# ------------------------------------------------------------
echo "=== Building GODAR ==="
mkdir -p "${GODAR_BUILD}"
cd "${GODAR_BUILD}"

cmake \
    -DCMAKE_Fortran_COMPILER=mpifort \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DCMAKE_PREFIX_PATH="${CORETRAN_INSTALL}" \
    -Dcoretran_DIR="${CORETRAN_INSTALL}/lib/cmake" \
    ..

make -j"$(nproc)"

echo ""
echo "============================================"
echo "  Build complete!"
echo "  GODAR binary: ${GODAR_SRC}/bin/godar"
echo "  Python venv:  ${VENV_DIR}"
echo "  Activate with: source ${VENV_DIR}/bin/activate"
echo "============================================"
