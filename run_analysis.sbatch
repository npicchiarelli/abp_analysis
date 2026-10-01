#!/usr/bin/env bash
# Run exec_analysis.jl on one run directory, copy the figures to processed_data
# and assemble the recap grids with make_recap.jl.
#
# Usage:  sbatch abp_analysis/run_analysis.sbatch <run_dir>
#   e.g.  sbatch -J ana_range abp_analysis/run_analysis.sbatch simulations/run02_range
#
# Outputs (for run_dir = .../run02_range):
#   <run_dir>/imgs/                                   figures from exec_analysis.jl
#   processed_data/run02_imgs/run02_range/imgs/       copy of those figures
#   processed_data/run02_imgs/recap/range/            recap_range_{polar,cluster,rdf,situa}.png
#
#SBATCH --job-name=abp_analysis
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.log
set -euo pipefail

ROOT=/storage/n.picchiarelli/abp_simulations
JULIA=/storage/n.picchiarelli/julia/usr/bin/julia
export JULIA_DEPOT_PATH=/storage/n.picchiarelli/.julia
export OPENBLAS_NUM_THREADS=1

RUN_DIR="$(realpath "${1:?usage: sbatch run_analysis.sbatch <run_dir>}")"
RUN_NAME="$(basename "$RUN_DIR")"          # run02_range
RUN_ID="${RUN_NAME%%_*}"                   # run02
METHOD="${RUN_NAME#*_}"                    # range
OUT="$ROOT/processed_data/${RUN_ID}_imgs"

echo "job $SLURM_JOB_ID on $(hostname)  cpus=$SLURM_CPUS_PER_TASK  $(date '+%F %T')"
echo "run dir: $RUN_DIR"

cd "$ROOT/abp_analysis"
"$JULIA" --project=. -t "$SLURM_CPUS_PER_TASK" exec_analysis.jl "$RUN_DIR"

mkdir -p "$OUT/$RUN_NAME/imgs"
cp "$RUN_DIR"/imgs/*.png "$OUT/$RUN_NAME/imgs/"

"$JULIA" --project=. make_recap.jl "$OUT/$RUN_NAME/imgs" "$OUT/recap/$METHOD" \
    --label "$RUN_ID, $METHOD" --prefix "$METHOD"

echo "done $(date '+%F %T')"
