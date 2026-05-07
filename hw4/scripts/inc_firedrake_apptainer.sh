# shellcheck shell=bash
# Source from hw4/scripts/*.sh — defines hw4_python3 for Firedrake runs on clusters.
#
# Set before running scripts:
#   export HW4_APPTAINER_SIF=/path/to/firedrake-ts.sif
#   export HW4_APPTAINER_BIND=/your/home:/your/home   # optional; default HOME:HOME
#
# Slurm + Open MPI in the image: isolated PLM is set unless already overridden.
# If you launch from `srun`/`sbatch` and see PMI / OPAL errors from MPI_Init, the
# container Open MPI is picking up SLURM_*; we default to unsetting those when
# SLURM_JOB_ID is set (override with HW4_APPTAINER_UNSET_SLURM=0 if your site needs PMI).

hw4_python3() {
  local sif="${HW4_APPTAINER_SIF:-}"
  if [[ -z "${sif}" && -n "${FIREDRAKE_TS_SIF:-}" ]]; then
    sif="${FIREDRAKE_TS_SIF}"
  fi

  if [[ -n "${sif}" ]]; then
    local bind="${HW4_APPTAINER_BIND:-${HOME}:${HOME}}"
    local runner="${HW4_APPTAINER_RUNNER:-}"
    if [[ -z "${runner}" ]]; then
      if command -v apptainer >/dev/null 2>&1; then
        runner=apptainer
      elif command -v singularity >/dev/null 2>&1; then
        runner=singularity
      else
        echo "ERROR: HW4_APPTAINER_SIF is set but neither apptainer nor singularity is in PATH." >&2
        return 1
      fi
    fi
    export OMPI_MCA_plm="${OMPI_MCA_plm:-isolated}"
    # Default: under Slurm, strip SLURM_* before `apptainer exec` so in-container
    # mpiexec does not attempt Slurm PMI (often missing / mismatched → MPI_Init failure).
    local unset_slurm="${HW4_APPTAINER_UNSET_SLURM:-}"
    if [[ -z "${unset_slurm}" && -n "${SLURM_JOB_ID:-}" ]]; then
      unset_slurm=1
    fi
    if [[ -z "${unset_slurm}" ]]; then
      unset_slurm=0
    fi
    if [[ "${unset_slurm}" == "1" ]]; then
      # shellcheck disable=SC2046
      unset $(printenv | grep '^SLURM_' | cut -d= -f1) || true
    fi
    # shellcheck disable=SC2086
    "${runner}" exec ${HW4_APPTAINER_EXTRA_ARGS:-} --bind "${bind}" "${sif}" python3 "$@"
  elif command -v mpiexec >/dev/null 2>&1; then
    mpiexec -n 1 python3 "$@"
  else
    python3 "$@"
  fi
}
