#!/bin/bash
set -euo pipefail

# ============================================================
# Build script for GODAR and its dependencies (coretran)
# on LUMI using amdflang via the mpifort wrapper.
#
# Run from any directory that contains src/ :
#   bash src/lumi/build.sh                  # incremental build
#   bash src/lumi/build.sh --clean-godar    # rebuild GODAR from scratch
#   bash src/lumi/build.sh --clean          # rebuild coretran and GODAR
# ============================================================

CLEAN_CORETRAN=0
CLEAN_GODAR=0
for _arg in "$@"; do
    case "${_arg}" in
        --clean)       CLEAN_CORETRAN=1; CLEAN_GODAR=1 ;;
        --clean-godar) CLEAN_GODAR=1 ;;
        -h|--help)     sed -n '4,15p' "$0"; exit 0 ;;
        *) echo "ERROR: unknown option ${_arg} (see --help)" >&2; exit 1 ;;
    esac
done
 
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
#    Only touch a file that needs it: rewriting it would change its
#    time stamp and force a recompilation.
for _f in "${TOPDIR}"/src/model/src/{godar,stepper,diagnostics_mpi,mpi_coms}.f90; do
    if grep -q 'use mpi_f08$' "$_f"; then
        sed -i 's/use mpi_f08$/use mpi/' "$_f"
    fi
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
# 3. Clean previous build artifacts (only when asked)
# ------------------------------------------------------------
if [ "${CLEAN_CORETRAN}" -eq 1 ]; then
    echo "=== Cleaning coretran build ==="
    rm -rf "${CORETRAN_SRC}/build" "${CORETRAN_INSTALL}"
fi
if [ "${CLEAN_GODAR}" -eq 1 ]; then
    echo "=== Cleaning GODAR build ==="
    rm -rf "${GODAR_BUILD}" "${GODAR_SRC}/bin" "${GODAR_SRC}/lib" "${GODAR_SRC}/include"
fi
echo ""
 
# ------------------------------------------------------------
# 4. Clone and build coretran (skipped if already installed)
# ------------------------------------------------------------
if [ -d "${CORETRAN_INSTALL}/lib/cmake" ]; then
    echo "=== coretran already installed in ${CORETRAN_INSTALL}: skipped ==="
    echo "    (use --clean to rebuild it)"
    echo ""
else
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
fi
 
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
    -DMPIEXEC_EXECUTABLE="$(command -v srun)" \
    ..
 
make -j"$(nproc)"
 
echo ""
echo "============================================"
echo "  Build complete!"
echo "  GODAR binary: ${GODAR_SRC}/bin/godar"
echo "  Python venv:  ${VENV_DIR}"
echo "  Activate with: source ${VENV_DIR}/bin/activate"
echo "============================================"