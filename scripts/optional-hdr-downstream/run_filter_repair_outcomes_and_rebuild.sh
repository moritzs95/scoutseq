#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${SCRIPTS_DIR}/.." && pwd)"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

resolve_path() {
  local path="$1"
  if [[ -z "${path}" ]]; then
    return
  fi
  if [[ "${path}" = /* ]]; then
    printf '%s\n' "${path}"
  else
    (cd "${CONFIG_DIR}" && python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${path}")
  fi
}

export_dro_target_config() {
  export SCOUT_DRO_LABEL="${DROTARGET:-}"
  export SCOUT_DRO_LEFT_ANCHOR="${DRO_LEFT_ANCHOR:-}"
  export SCOUT_DRO_RIGHT_ANCHOR="${DRO_RIGHT_ANCHOR:-}"
  export SCOUT_DRO_HDR_ANCHOR="${DRO_HDR_ANCHOR:-}"
  export SCOUT_DRO_WT_AMPLICON="${DRO_WT_AMPLICON:-}"
  export SCOUT_DRO_WT_AMPLICON_AFTER_CUTSITE="${DRO_WT_AMPLICON_AFTER_CUTSITE:-}"
  export SCOUT_DRO_HDRBC_LENGTH="${HDRBCLENGTH:-}"
  export SCOUT_DRO_GS_PRIMER="${DRO_GS_PRIMER:-}"
  export SCOUT_DRO_CDNA_LENGTH="${DRO_CDNA_LENGTH:-}"
  export SCOUT_DRO_REFERENCE_LENGTH="${DRO_REFERENCE_LENGTH:-}"
  export SCOUT_DRO_INSERTION_TYPE="${DRO_INSERTION_TYPE:-}"
  export SCOUT_DRO_CINS_HDRBC="${DRO_CINS_HDRBC:-}"
  export SCOUT_DRO_THRESHOLD_I90PLUS="${DRO_THRESHOLD_I90PLUS:-}"
  export SCOUT_DRO_THRESHOLD_D90PLUS="${DRO_THRESHOLD_D90PLUS:-}"
  export SCOUT_DRO_THRESHOLD_S1PLUS="${DRO_THRESHOLD_S1PLUS:-}"
}

CONFIG_PATH="${1:-}"
OUTCOMES_LIST="${2:-}"
REPAIR_TABLE_INPUT="${3:-}"
OUTPUT_DIR_INPUT="${4:-}"

[[ -n "${CONFIG_PATH}" && -f "${CONFIG_PATH}" ]] || die "Usage: $0 <config-file> <outcomes-list> [repair-table] [output-dir]"
[[ -n "${OUTCOMES_LIST}" && -f "${OUTCOMES_LIST}" ]] || die "Outcome list '${OUTCOMES_LIST}' does not exist."

CONFIG_PATH="$(cd "$(dirname "${CONFIG_PATH}")" && pwd)/$(basename "${CONFIG_PATH}")"
CONFIG_DIR="$(dirname "${CONFIG_PATH}")"
source "${CONFIG_PATH}"

export SCOUT_LIBTYPE="${LIBTYPE:-10X}"
export_dro_target_config

if [[ -n "${OUTPUTDIR:-}" ]]; then
  RUN_DIR="$(resolve_path "${OUTPUTDIR}")"
else
  [[ -n "${SAMPLENAME:-}" ]] || die "SAMPLENAME is required in the config."
  RUN_DIR="${ROOT_DIR}/${SAMPLENAME}"
fi

if [[ -n "${REPAIR_TABLE_INPUT}" ]]; then
  REPAIR_TABLE="$(resolve_path "${REPAIR_TABLE_INPUT}")"
else
  REPAIR_TABLE="${RUN_DIR}/$(basename "${RUN_DIR}")_crispresso_out"
  if [[ -d "${REPAIR_TABLE}" ]]; then
    die "Automatic repair-table detection is not supported for this layout. Pass the RepairOutcomeTable.csv path explicitly."
  fi
  REPAIR_TABLE="${PWD}/RepairOutcomeTable.csv"
fi

[[ -f "${REPAIR_TABLE}" ]] || die "RepairOutcomeTable.csv '${REPAIR_TABLE}' does not exist."

if [[ -n "${OUTPUT_DIR_INPUT}" ]]; then
  OUTPUT_DIR="$(resolve_path "${OUTPUT_DIR_INPUT}")"
else
  OUTPUT_DIR="$(cd "$(dirname "${REPAIR_TABLE}")" && pwd)/filtered_repair_outcomes"
fi

CBC_DIR="$(resolve_path "${CBCPATH:-}")"
[[ -d "${CBC_DIR}" ]] || die "CBCPATH '${CBC_DIR}' does not exist."

TRANSLOCATION_SUMMARY=""
REPAIR_DIR="$(cd "$(dirname "${REPAIR_TABLE}")" && pwd)"
if [[ -f "${REPAIR_DIR}/../TranslocationSelection.csv" ]]; then
  TRANSLOCATION_SUMMARY="${REPAIR_DIR}/../TranslocationSelection.csv"
elif [[ -f "${REPAIR_DIR}/../translocation_coordinates_filtered_output.csv" ]]; then
  TRANSLOCATION_SUMMARY="${REPAIR_DIR}/../translocation_coordinates_filtered_output.csv"
fi

python3 "${SCRIPT_DIR}/filter_repair_outcomes_and_rebuild.py" \
  --repair-table "${REPAIR_TABLE}" \
  --outcomes-list "$(resolve_path "${OUTCOMES_LIST}")" \
  --cbc-dir "${CBC_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --libtype "${LIBTYPE:-10X}" \
  --target-chromosome "${CHRTARGET:-}" \
  --translocations "${TRANSLOCATION_SUMMARY}"
