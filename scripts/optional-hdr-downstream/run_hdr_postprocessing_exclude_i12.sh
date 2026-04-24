#!/usr/bin/env bash

# Standalone downstream HDR postprocessing runner that excludes CRISPResso reads
# containing I12 insertions before running the standard downstream HDR scripts.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${SCRIPTS_DIR}/.." && pwd)"
HDR_DIR="${SCRIPTS_DIR}/optional-hdr-downstream"
PARSE_DIR="${SCRIPTS_DIR}/optional-parse"
GENOME_DIR_DEFAULT="${ROOT_DIR}/genome"

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

log() {
  printf '%s [scOUT-i12] %s\n' "$(timestamp)" "$*"
}

die() {
  printf '%s [scOUT-i12] ERROR: %s\n' "$(timestamp)" "$*" >&2
  exit 1
}

require_var() {
  local var_name="$1"
  if [[ -z "${!var_name:-}" ]]; then
    die "Required config variable '${var_name}' is missing."
  fi
}

resolve_path() {
  local value="$1"
  if [[ -z "${value}" ]]; then
    return 0
  fi
  if [[ "${value}" = /* ]]; then
    printf '%s\n' "${value}"
  else
    printf '%s\n' "${CONFIG_DIR}/${value}"
  fi
}

export_dro_target_config() {
  export SCOUT_DRO_LABEL="${DROTARGET:-}"
  require_var "DRO_LEFT_ANCHOR"
  require_var "DRO_RIGHT_ANCHOR"
  require_var "DRO_WT_AMPLICON"
  require_var "DRO_GS_PRIMER"
  require_var "DRO_CDNA_LENGTH"
  require_var "DRO_REFERENCE_LENGTH"

  export SCOUT_DRO_LEFT_ANCHOR="${DRO_LEFT_ANCHOR}"
  export SCOUT_DRO_RIGHT_ANCHOR="${DRO_RIGHT_ANCHOR}"
  export SCOUT_DRO_HDR_ANCHOR="${DRO_HDR_ANCHOR:-}"
  export SCOUT_DRO_WT_AMPLICON="${DRO_WT_AMPLICON}"
  export SCOUT_DRO_WT_AMPLICON_AFTER_CUTSITE="${DRO_WT_AMPLICON_AFTER_CUTSITE:-${DRO_WT_AMPLICON}}"
  export SCOUT_DRO_HDRBC_LENGTH="${HDRBCLENGTH}"
  export SCOUT_DRO_GS_PRIMER="${DRO_GS_PRIMER}"
  export SCOUT_DRO_CDNA_LENGTH="${DRO_CDNA_LENGTH}"
  export SCOUT_DRO_REFERENCE_LENGTH="${DRO_REFERENCE_LENGTH}"
  export SCOUT_DRO_THRESHOLD_I90PLUS="${DRO_THRESHOLD_I90PLUS:-30}"
  export SCOUT_DRO_THRESHOLD_D90PLUS="${DRO_THRESHOLD_D90PLUS:-30}"
  export SCOUT_DRO_THRESHOLD_S1PLUS="${DRO_THRESHOLD_S1PLUS:-50}"
  export SCOUT_DRO_INSERTION_TYPE="${DRO_INSERTION_TYPE:-substitution}"

  if [[ "${SCOUT_DRO_INSERTION_TYPE}" != "rINS" ]]; then
    require_var "DRO_HDR_ANCHOR"
  fi

  if [[ "${SCOUT_DRO_INSERTION_TYPE}" == "cINS" ]]; then
    require_var "DRO_CINS_HDRBC"
    export SCOUT_DRO_CINS_HDRBC="${DRO_CINS_HDRBC}"
  else
    unset SCOUT_DRO_CINS_HDRBC || true
  fi
}

main() {
  local config_path="${1:-}"
  [[ -n "${config_path}" ]] || die "Usage: $0 <config-file>"
  [[ -f "${config_path}" ]] || die "Config file '${config_path}' does not exist."

  CONFIG_PATH="$(cd "$(dirname "${config_path}")" && pwd)/$(basename "${config_path}")"
  CONFIG_DIR="$(dirname "${CONFIG_PATH}")"

  # shellcheck disable=SC1090
  source "${CONFIG_PATH}"

  require_var "SAMPLENAME"
  require_var "LIBTYPE"
  require_var "FILTERHDRREADSCONFIG"
  require_var "HDRBCLENGTH"
  require_var "DROREADLENGTH"
  require_var "CBCPATH"
  require_var "SPECIES"
  require_var "CHRTARGET"
  require_var "STARTTARGET"
  require_var "ENDTARGET"
  require_var "I12_FILTER_LEFT_ANCHOR"
  require_var "I12_FILTER_RIGHT_ANCHOR"
  require_var "I12_FILTER_BP_BETWEEN_ANCHORS"

  CBCPATH="$(resolve_path "${CBCPATH}")"
  GENOME_DIR="$(resolve_path "${GENOME_DIR:-${GENOME_DIR_DEFAULT}}")"
  AAV_GENOME_DIR="$(resolve_path "${AAV_GENOME_DIR:-${GENOME_DIR}}")"

  local run_dir
  if [[ -n "${OUTPUTDIR:-}" ]]; then
    run_dir="$(resolve_path "${OUTPUTDIR}")"
  else
    run_dir="${ROOT_DIR}/${SAMPLENAME}"
  fi

  local fastq_corrected_dir="${run_dir}/fastq_corrected_bc"
  local crispresso_input_r2="${fastq_corrected_dir}/${SAMPLENAME}.extracted.noBC.R2.fastq.gz"
  local hdrbc_r2="${fastq_corrected_dir}/${SAMPLENAME}.extracted.HDRBC.R2.fastq.gz"
  if [[ -z "${FILTERHDRREADSCONFIG:-}" ]]; then
    crispresso_input_r2="${fastq_corrected_dir}/${SAMPLENAME}.extracted.R2.fastq.gz"
  fi

  local crispresso_dir_name
  crispresso_dir_name="$(basename "${crispresso_input_r2}" .fastq.gz)"
  local crispresso_result_dir="${run_dir}/${SAMPLENAME}_crispresso_out/CRISPResso_on_${crispresso_dir_name}"

  [[ -d "${crispresso_result_dir}" ]] || die "CRISPResso result directory '${crispresso_result_dir}' does not exist."
  [[ -f "${crispresso_result_dir}/CRISPResso_output.fastq.gz" ]] || die "CRISPResso output FASTQ is missing in '${crispresso_result_dir}'."

  export SCOUT_LIBTYPE="${LIBTYPE:-}"
  export_dro_target_config

  cd "${crispresso_result_dir}"
  log "Filtering CRISPResso output to exclude reads with I12 insertions."

  local original_crispresso="CRISPResso_output.fastq.gz"
  local backup_crispresso="CRISPResso_output.with_i12.fastq.gz"
  local filtered_crispresso="CRISPResso_output.no_i12.fastq.gz"
  local i12_max_mismatches="${I12_FILTER_MAX_MISMATCHES:-1}"

  rm -f "${filtered_crispresso}"
  python "${HDR_DIR}/filter_crispresso_exclude_i12.py" \
    "${original_crispresso}" \
    "${filtered_crispresso}" \
    "${I12_FILTER_LEFT_ANCHOR}" \
    "${I12_FILTER_RIGHT_ANCHOR}" \
    "${I12_FILTER_BP_BETWEEN_ANCHORS}" \
    "${i12_max_mismatches}"
  mv "${original_crispresso}" "${backup_crispresso}"
  trap 'if [[ -f "'"${backup_crispresso}"'" ]]; then mv -f "'"${backup_crispresso}"'" "'"${original_crispresso}"'"; fi' EXIT
  cp "${filtered_crispresso}" "${original_crispresso}"

  log "Running standard downstream HDR postprocessing on the filtered CRISPResso output."

  python "${HDR_DIR}/crispresso_fastq_to_table.py" \
    "${hdrbc_r2}" \
    "${HDRBCLENGTH}" \
    "${LIBTYPE}" \
    "${DROTARGET}" \
    "${DROREADLENGTH}" \
    "${PARSEKIT:-}"

  SCOUT_GENOME_DIR="${GENOME_DIR}" bash "${HDR_DIR}/offtargets_bwa_mem.sh" \
    "${CHRTARGET}" \
    "${STARTTARGET}" \
    "${ENDTARGET}" \
    potential_offtargets.csv \
    "${SPECIES}" \
    "${LIBTYPE}"

  SCOUT_GENOME_DIR="${GENOME_DIR}" bash "${HDR_DIR}/translocations_bwa_mem.sh" \
    "${CHRTARGET}" \
    "${STARTTARGET}" \
    "${ENDTARGET}" \
    potential_translocations.csv \
    "${SPECIES}" \
    "${LIBTYPE}"

  if [[ "${LIBTYPE}" == "PARSE" ]]; then
    python "${PARSE_DIR}/Parse_translocation_extractor.py" "${CBCPATH}" translocation_coordinates_filtered.csv
  else
    python "${HDR_DIR}/10X_translocation_extractor.py"
  fi

  if [[ -f unaligned_translocations.csv ]] && [[ "$(wc -l < unaligned_translocations.csv)" -gt 1 ]]; then
    SCOUT_GENOME_DIR="${AAV_GENOME_DIR}" bash "${HDR_DIR}/unaligned_translocation_bwa_mem.sh" unaligned_translocations.csv "${LIBTYPE}"
    if [[ "${LIBTYPE}" == "PARSE" ]]; then
      python "${PARSE_DIR}/Parse_translocation_extractor.py" "${CBCPATH}" AAV_translocation_coordinates.csv
    else
      python "${HDR_DIR}/10X_translocation_extractor_AAV.py" AAV_translocation_coordinates.csv
    fi
  else
    log "No non-empty unaligned translocation file detected; skipping AAV remapping."
  fi

  mkdir -p editing_outcomes
  cd editing_outcomes

  if [[ "${LIBTYPE}" == "PARSE" ]]; then
    python -u "${PARSE_DIR}/Parse_scOUT_seqDRO_extractor.py" \
      "${DROTARGET}" \
      "${DROREADLENGTH}" \
      "${CBCPATH}" \
      > EditingOutcomeAssignment.log 2>&1
  else
    python -u "${HDR_DIR}/10X_scOUT_seqDRO_extractor.py" \
      "${DROTARGET}" \
      "${DROREADLENGTH}" \
      "${CBCPATH}" \
      > EditingOutcomeAssignment.log 2>&1
  fi

  local translocation_summary=""
  if [[ -f ../TranslocationSelection.csv ]]; then
    translocation_summary="../TranslocationSelection.csv"
  elif [[ -f ../translocation_coordinates_filtered_output.csv ]]; then
    translocation_summary="../translocation_coordinates_filtered_output.csv"
  fi

  python "${SCRIPTS_DIR}/required-core/integrate_translocations_and_filter_outcomes.py" \
    --editing-outcomes EditingOutcomesByCellBarcode.xlsx \
    --output-dir "${PWD}" \
    --target-chromosome "${CHRTARGET:-}" \
    --translocations "${translocation_summary}"

  mv -f "${backup_crispresso}" "${original_crispresso}"
  trap - EXIT
  log "Downstream HDR postprocessing with I12 exclusion completed."
}

main "$@"
