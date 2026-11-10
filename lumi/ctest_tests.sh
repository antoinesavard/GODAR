#!/bin/bash
#SBATCH --job-name=godar-ctest
#SBATCH --account=project_462001545
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2              # largest run in cases.json: mpi2
#SBATCH --cpus-per-task=4       # largest run in cases.json: omp4
#SBATCH --mem=16G               # "several GB" per case
#SBATCH --output=src/lumi/%j.out

set -euo pipefail

# directory sbatch was run from: the parent of src/
export TOPDIR="${SLURM_SUBMIT_DIR}"
BUILD="${TOPDIR}/src/build"

echo "TOPDIR = ${TOPDIR}"
if [ ! -d "${BUILD}" ]; then
    echo "No build directory at ${BUILD}: submit from the parent of src/" >&2
    exit 1
fi

source "${TOPDIR}/src/lumi/env.sh"
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-4}

cd "${BUILD}"
ctest --output-on-failure -V
