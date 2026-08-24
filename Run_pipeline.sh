#!/usr/bin/env bash
# ===========================================================================
# Run_pipeline.sh — robust entry point for the UV-RIP Snakemake workflow.
#
# Handles the two things that reliably break Snakemake on Quest:
#   1. conda not being a shell function in non-interactive shells
#      (sources conda.sh from $CONDA_EXE, with fallbacks).
#   2. system modules (e.g. deeptools) injecting the wrong Python
#      (runs `module purge` after activating the conda env).
#
# Usage:
#   bash Run_pipeline.sh [--dry-run] [--slurm] [--use-conda] [--unlock] \
#                        [-- <extra snakemake args>]
#
#   bash Run_pipeline.sh --dry-run            # local dry run
#   bash Run_pipeline.sh --slurm --dry-run    # SLURM dry run (uses profile)
#   bash Run_pipeline.sh --slurm              # full SLURM run
#   bash Run_pipeline.sh                      # full local run
# ===========================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CONDA_ENV="${UVRIP_ENV:-uvrip}"
CORES="${CORES:-20}"

DRY=""; SLURM=0; USECONDA=""; UNLOCK=0; PASSTHRU=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run|-n) DRY="-n" ;;
    --slurm)      SLURM=1 ;;
    --use-conda)  USECONDA="--use-conda --conda-frontend mamba" ;;
    --unlock)     UNLOCK=1 ;;
    --)           shift; PASSTHRU=("$@"); break ;;
    *)            PASSTHRU+=("$1") ;;
  esac
  shift
done

# ---- Portable conda activation -------------------------------------------
CONDA_SH=""
if [[ -n "${CONDA_EXE:-}" ]]; then
  CONDA_SH="$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
fi
for c in "$CONDA_SH" \
         "$HOME/miniconda3/etc/profile.d/conda.sh" \
         "$HOME/miniforge3/etc/profile.d/conda.sh" \
         "$HOME/.conda/etc/profile.d/conda.sh" \
         "/hpc/software/mamba/23.1.0/etc/profile.d/conda.sh"; do
  if [[ -n "$c" && -f "$c" ]]; then CONDA_SH="$c"; break; fi
done
[[ -f "$CONDA_SH" ]] || { echo "ERROR: cannot find conda.sh; set CONDA_EXE." >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$CONDA_ENV"
echo "conda env active: $CONDA_ENV  (from $CONDA_SH)"

# Clear any inherited modules so system Python/tools don't shadow the env.
module purge 2>/dev/null || true

# Help --use-conda subshells find conda (only matters if --use-conda is passed).
export PATH="$(dirname "$(dirname "$CONDA_SH")")/../condabin:$PATH"

command -v snakemake >/dev/null || { echo "ERROR: snakemake not on PATH in env $CONDA_ENV." >&2; exit 1; }
echo "snakemake: $(snakemake --version)"

[[ "$UNLOCK" -eq 1 ]] && snakemake --unlock

if [[ "$SLURM" -eq 1 ]]; then
  snakemake --profile profiles/slurm $USECONDA $DRY --rerun-incomplete "${PASSTHRU[@]}"
else
  snakemake --cores "$CORES" $USECONDA $DRY --rerun-incomplete "${PASSTHRU[@]}"
fi
