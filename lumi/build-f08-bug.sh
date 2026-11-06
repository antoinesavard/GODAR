#!/bin/bash
set -euo pipefail

# ============================================================
# Build GODAR with mpi_f08 to reproduce the CFI descriptor bug.
#
# Patches `use mpi` -> `use mpi_f08` in the 4 source files,
# rebuilds, and submits a quick 1-config test job.
#
# Usage:  bash src/lumi/build-f08-bug.sh
# Revert: bash src/lumi/build.sh   (rebuilds with use mpi)
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOPDIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
GODAR_SRC="${TOPDIR}/src"

SOURCES=(
    "${GODAR_SRC}/model/src/godar.f90"
    "${GODAR_SRC}/model/src/stepper.f90"
    "${GODAR_SRC}/model/src/diagnostics_mpi.f90"
    "${GODAR_SRC}/model/src/mpi_coms.f90"
)

# ============================================================
# 1. Patch: use mpi  ->  use mpi_f08
#    (only lines ending with exactly "use mpi", not "use mpi_var" etc.)
# ============================================================
echo "=== Patching sources: use mpi -> use mpi_f08 ==="
for f in "${SOURCES[@]}"; do
    before=$(grep -c 'use mpi$' "$f" || true)
    sed -i 's/use mpi$/use mpi_f08/' "$f"
    after=$(grep -c 'use mpi_f08' "$f" || true)
    echo "  $(basename "$f"): ${before} -> ${after} occurrence(s)"
done
echo ""

# ============================================================
# 2. Build (reuses the normal build script's cmake logic)
# ============================================================
echo "=== Building GODAR (mpi_f08 variant) ==="
source "${SCRIPT_DIR}/build.sh"
echo ""

# ============================================================
# 3. Submit a quick 1-config test (1 node, 2 ranks, 7 threads)
# ============================================================
WORK="${GODAR_SRC}/tests/work/f08_bug"
mkdir -p "${WORK}"

cat > "${WORK}/job.sh" << JOBEOF
#!/bin/bash
#SBATCH --job-name=godar-f08-bug
#SBATCH --account=project_462000125
#SBATCH --partition=standard-g
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:15:00
#SBATCH --output=${WORK}/slurm-%j.out
#SBATCH --error=${WORK}/slurm-%j.err

export TOPDIR="${TOPDIR}"
source "${SCRIPT_DIR}/env.sh"

cd "${GODAR_SRC}"
python3 tests/run_bench.py \\
    --particles 20000 \\
    --nodes 1 \\
    --ranks-per-node 2 \\
    --threads 7 \\
    --repeat 1 \\
    --work "${WORK}"
JOBEOF

echo "=== Submitting f08-bug test job ==="
sbatch "${WORK}/job.sh"
echo ""
echo "Monitor with:  tail -f ${WORK}/slurm-*.out"
echo "Check log:     cat ${WORK}/ARCTIC_n20000/log.candidate"
