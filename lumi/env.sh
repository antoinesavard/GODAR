# ============================================================
# LUMI environment setup — sourced by job scripts and by
# run_benchmark.sh.  Also usable interactively:
#
#   source src/lumi/env.sh
#
# Expects TOPDIR to be set (the repo root, parent of src/).
# If not set, it is derived from this script's location.
# ============================================================

if [ -z "${TOPDIR:-}" ]; then
    TOPDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
fi

module use /appl/local/containers/test-modules
ml rocm-afar-cray-mpich/24.3.0
ml netcdf_fortran
ml netcdf_c
ml cray-python

source "${TOPDIR}/venv/bin/activate"

# Increase OpenMP thread stack size for recursive KD-tree searches
export OMP_STACKSIZE=64M
