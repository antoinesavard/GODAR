#!/bin/bash
#SBATCH --job-name=godar-bench
#SBATCH --account=project_462001545
#SBATCH --partition=standard    # CPU nodes, 128 cores
#SBATCH --nodes=4               # largest number of nodes of the benchmark
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=src/lumi/%j.out

# ============================================================
# ARCTIC benchmark of GODAR, run by ctest. After building GODAR,
# submit from the parent of src/:
#   sbatch src/lumi/ctest_bench.sh
#
# ctest builds the ARCTIC case (test arctic_case, kept for the next
# runs), then times godar on 1, 2, 4 nodes x 1, 2, 4 MPI processes
# per node x 16 to 128 threads (test arctic, tests/CMakeLists.txt).
# The configurations that do not fit the allocation are skipped. The
# number of floes is set when configuring: -DGODAR_BENCH_PARTICLES
# (default 20000). The results are also appended to
# src/tests/work/bench/bench.csv.
# ============================================================

set -euo pipefail

if [ -z "${SLURM_JOB_ID:-}" ]; then
    echo "Submit this with: sbatch $0" >&2
    exit 1
fi

# directory sbatch was run from: the parent of src/
export TOPDIR="${SLURM_SUBMIT_DIR}"
source "${TOPDIR}/src/lumi/env.sh"

cd "${TOPDIR}/src/build"
ctest -C bench -R arctic -V
