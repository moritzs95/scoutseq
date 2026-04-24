#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${SCRIPTS_DIR}/.." && pwd)"
PARSE_DIR="${SCRIPTS_DIR}/optional-parse"
HDR_DIR="${SCRIPTS_DIR}/optional-hdr-downstream"
CORE_DIR="${SCRIPTS_DIR}/required-core"

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
  export SCOUT_DRO_WT_AMPLICON_AFTER_CUTSITE="${DRO_WT_AMPLICON_AFTER_CUTSITE:-${DRO_WT_AMPLICON:-}}"
  export SCOUT_DRO_HDRBC_LENGTH="${HDRBCLENGTH:-}"
  export SCOUT_DRO_GS_PRIMER="${DRO_GS_PRIMER:-}"
  export SCOUT_DRO_CDNA_LENGTH="${DRO_CDNA_LENGTH:-}"
  export SCOUT_DRO_REFERENCE_LENGTH="${DRO_REFERENCE_LENGTH:-}"
  export SCOUT_DRO_INSERTION_TYPE="${DRO_INSERTION_TYPE:-}"
  export SCOUT_DRO_CINS_HDRBC="${DRO_CINS_HDRBC:-}"
  export SCOUT_DRO_THRESHOLD_I90PLUS="${DRO_THRESHOLD_I90PLUS:-30}"
  export SCOUT_DRO_THRESHOLD_D90PLUS="${DRO_THRESHOLD_D90PLUS:-30}"
  export SCOUT_DRO_THRESHOLD_S1PLUS="${DRO_THRESHOLD_S1PLUS:-50}"
}

CONFIG_PATH="${1:-}"
MODS_TABLE_INPUT="${2:-}"
OUTPUT_ROOT_INPUT="${3:-}"

[[ -n "${CONFIG_PATH}" && -f "${CONFIG_PATH}" ]] || die "Usage: $0 <config-file> [bc_umi_mod_seq.csv] [output-root]"

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

FASTQ_CORRECTED_DIR="${RUN_DIR}/fastq_corrected_bc"
CRISPRESSO_INPUT_R2="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.R2.fastq.gz"
if [[ -n "${FILTERHDRREADSCONFIG:-}" ]]; then
  CRISPRESSO_INPUT_R2="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.noBC.R2.fastq.gz"
fi

CRISPRESSO_DIR_NAME="$(basename "${CRISPRESSO_INPUT_R2}" .fastq.gz)"
CRISPRESSO_RESULT_DIR="${RUN_DIR}/${SAMPLENAME}_crispresso_out/CRISPResso_on_${CRISPRESSO_DIR_NAME}"

if [[ -n "${MODS_TABLE_INPUT}" ]]; then
  MODS_TABLE="$(resolve_path "${MODS_TABLE_INPUT}")"
else
  MODS_TABLE="${CRISPRESSO_RESULT_DIR}/bc_umi_mod_seq.csv"
fi
[[ -f "${MODS_TABLE}" ]] || die "bc_umi_mod_seq.csv '${MODS_TABLE}' does not exist."

if [[ -n "${OUTPUT_ROOT_INPUT}" ]]; then
  OUTPUT_ROOT="$(resolve_path "${OUTPUT_ROOT_INPUT}")"
else
  OUTPUT_ROOT="${CRISPRESSO_RESULT_DIR}/non_hdr_only"
fi

WORK_ROOT="${OUTPUT_ROOT}"
EDITING_OUTCOMES_DIR="${WORK_ROOT}/editing_outcomes"
HDR_MODS_PATH="${WORK_ROOT}/Hdr_mods.csv"
mkdir -p "${EDITING_OUTCOMES_DIR}"

CBCPATH_RESOLVED="$(resolve_path "${CBCPATH:-}")"
[[ -d "${CBCPATH_RESOLVED}" ]] || die "CBCPATH '${CBCPATH_RESOLVED}' does not exist."

TRANSLOCATION_SUMMARY=""
if [[ -f "${CRISPRESSO_RESULT_DIR}/TranslocationSelection.csv" ]]; then
  TRANSLOCATION_SUMMARY="${CRISPRESSO_RESULT_DIR}/TranslocationSelection.csv"
elif [[ -f "${CRISPRESSO_RESULT_DIR}/translocation_coordinates_filtered_output.csv" ]]; then
  TRANSLOCATION_SUMMARY="${CRISPRESSO_RESULT_DIR}/translocation_coordinates_filtered_output.csv"
fi

python3 "${SCRIPT_DIR}/build_non_hdr_mods_table.py" \
  --mods-table "${MODS_TABLE}" \
  --output "${HDR_MODS_PATH}"

cd "${EDITING_OUTCOMES_DIR}"

if [[ "${LIBTYPE:-}" == "PARSE" ]]; then
  python -u "${PARSE_DIR}/Parse_scOUT_seqDRO_extractor.py" \
    "${DROTARGET}" \
    "${DROREADLENGTH}" \
    "${CBCPATH_RESOLVED}" \
    > EditingOutcomeAssignment.log 2>&1
else
  python -u "${HDR_DIR}/10X_scOUT_seqDRO_extractor.py" \
    "${DROTARGET}" \
    "${DROREADLENGTH}" \
    "${CBCPATH_RESOLVED}" \
    > EditingOutcomeAssignment.log 2>&1
fi

python3 "${SCRIPT_DIR}/filter_incorrect_hdr_from_editing_outcomes.py" \
  --input "${EDITING_OUTCOMES_DIR}/EditingOutcomesByCellBarcode.xlsx" \
  --output "${EDITING_OUTCOMES_DIR}/EditingOutcomesByCellBarcode.xlsx"

python3 "${CORE_DIR}/integrate_translocations_and_filter_outcomes.py" \
  --editing-outcomes EditingOutcomesByCellBarcode.xlsx \
  --output-dir "${EDITING_OUTCOMES_DIR}" \
  --target-chromosome "${CHRTARGET:-}" \
  --translocations "${TRANSLOCATION_SUMMARY}"

printf 'Non-HDR-only outcomes written to %s\n' "${EDITING_OUTCOMES_DIR}"
