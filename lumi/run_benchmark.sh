#!/bin/bash
set -euo pipefail

# ============================================================
# GODAR benchmark launcher for LUMI standard-g
#
# Run from any directory that contains src/ :
#   bash src/lumi/run_benchmark.sh
#
# Step 1: pre-build the ARCTIC case on the login node
#         (needs internet + numpy/scipy/numba/Pillow)
# Step 2: generate and submit a Slurm job script via
#         run_bench.py --slurm-script
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOPDIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [ ! -d "${TOPDIR}/src" ]; then
    echo "ERROR: cannot find src/ under ${TOPDIR}" >&2
    exit 1
fi

GODAR_SRC="${TOPDIR}/src"
ACCOUNT="project_462000125"
PARTITION="standard-g"
WALLTIME="2:00:00"
NODES="1,2,4,8"
RANKS_PER_NODE="1,2"
THREADS="7,14,28,56"
PARTICLES=20000
JOB_SCRIPT="${TOPDIR}/bench_job.sh"

# ------------------------------------------------------------
# 1. Load modules and activate venv
# ------------------------------------------------------------
echo "=== Loading environment ==="
export TOPDIR
source "${SCRIPT_DIR}/env.sh"

echo "Python: $(which python3) — $(python3 --version)"
echo "ncdump: $(which ncdump)"
echo ""

# ------------------------------------------------------------
# 2. Install Python dependencies for ARCTIC case generation
# ------------------------------------------------------------
echo "=== Installing Python dependencies for ARCTIC case ==="
pip install --quiet numpy scipy numba Pillow matplotlib
echo ""

# ------------------------------------------------------------
# 3. Pre-build the ARCTIC case on the login node
#    (compute nodes may lack internet and these tools)
# ------------------------------------------------------------
echo "=== Pre-building ARCTIC case with ${PARTICLES} floes ==="
cd "${GODAR_SRC}"
python3 tests/bench/make_arctic.py \
    --n "${PARTICLES}" \
    --out tests/work/bench/cases/ARCTIC_${PARTICLES}
echo ""

# ------------------------------------------------------------
# 4. Generate the Slurm job script via run_bench.py
# ------------------------------------------------------------
echo "=== Generating Slurm job script ==="
python3 tests/run_bench.py \
    --particles "${PARTICLES}" \
    --nodes "${NODES}" \
    --ranks-per-node "${RANKS_PER_NODE}" \
    --threads "${THREADS}" \
    --slurm-script "${JOB_SCRIPT}" \
    --account "${ACCOUNT}" \
    --partition "${PARTITION}" \
    --time "${WALLTIME}" \
    --setup "${SCRIPT_DIR}/env.sh"
echo ""

# ------------------------------------------------------------
# 5. Submit the job
# ------------------------------------------------------------
echo "=== Submitting job ==="
echo "Job script: ${JOB_SCRIPT}"
cat "${JOB_SCRIPT}"
echo ""
echo "---"
sbatch "${JOB_SCRIPT}"
