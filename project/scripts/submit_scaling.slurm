#!/bin/bash
# Optional Slurm driver: edit modules / paths for your site. Primary workflow is local PETSc
# (see docs/RUN_LATER_CLUSTER_APPTAINER.md and README.md).
#SBATCH --job-name=mc-scaling
#SBATCH --output=output/slurm-%j.out
#SBATCH --error=output/slurm-%j.err
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# Edit for your cluster environment:
# module load python
# module load petsc
# module load openmpi

set -euo pipefail

cd "$(dirname "$0")/.."

if [ ! -d ".venv" ]; then
  python -m venv .venv
fi
source .venv/bin/activate
pip install -r requirements.txt

python scripts/run_scaling.py \
  --mode strong \
  --ranks 1 2 4 8 \
  --base-paths 400000 \
  --option european \
  --out output/strong_scaling.csv

python scripts/plot_scaling.py \
  --infile output/strong_scaling.csv \
  --out-prefix output/strong_scaling
