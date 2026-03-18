#!/usr/bin/env bash

# Main entrypoint for the scOUTseq analysis workflow. This script stitches together FASTQ preprocessing, barcode extraction and salvaging, CRISPResso alignment and editing outcome categorization, optional HDR barcode handling, translocation/large deletion detection and assignment, and the final editing outcome reports.

set -euo pipefail

# Resolve the workspace layout once so the pipeline can run from any directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${SCRIPTS_DIR}/.." && pwd)"
GENOME_DIR_DEFAULT="${ROOT_DIR}/genome"
PARSE_DIR="${SCRIPTS_DIR}/optional-parse"
HDR_DIR="${SCRIPTS_DIR}/optional-hdr-downstream"
LOG_FILE=""

timestamp() {
  date '+%Y-%m-%d %H:%M:%S'
}

log() {
  printf '%s [scOUT] %s\n' "$(timestamp)" "$*"
}

die() {
  printf '%s [scOUT] ERROR: %s\n' "$(timestamp)" "$*" >&2
  exit 1
}

setup_logging() {
  local log_dir
  log_dir="${RUN_DIR}/logs"
  mkdir -p "${log_dir}"

  LOG_FILE="${log_dir}/${SAMPLENAME}_$(date '+%Y%m%d_%H%M%S').log"
  touch "${LOG_FILE}"

  # Mirror all subsequent stdout/stderr to a per-run log file.
  exec > >(tee -a "${LOG_FILE}") 2>&1
  log "Writing run log to ${LOG_FILE}."
}

load_conda() {
  if ! command -v conda >/dev/null 2>&1; then
    die "conda is required but was not found in PATH."
  fi

  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
}

activate_conda_if_available() {
  local env_name="$1"

  if conda env list | awk '{print $1}' | grep -Fxq "${env_name}"; then
    log "Activating conda environment '${env_name}'."
    conda activate "${env_name}"
  else
    log "Conda environment '${env_name}' is not available. Continuing with current environment."
  fi
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

run_seq_filter() {
  # Optionally drop reads with known sequence motifs before barcode handling.
  if [[ -n "${SEQFILTER:-}" ]]; then
    local seq_filter_path
    seq_filter_path="$(resolve_path "${SEQFILTER}")"
    log "Filtering reads with patterns from ${seq_filter_path}."

    seqtk mergepe "${RAW_R1}" "${RAW_R2}" \
      | python "${SCRIPT_DIR}/filter_fastq_with_patterns.py" "${seq_filter_path}" \
      | tee >(seqtk seq -1 - | pigz -c > "${CLEAN_R1}") \
      | seqtk seq -2 - | pigz -c > "${CLEAN_R2}"
  else
    ln -sfn "${RAW_R1}" "${CLEAN_R1}"
    ln -sfn "${RAW_R2}" "${CLEAN_R2}"
  fi
}

prepare_hdr_filter_config() {
  if [[ -z "${FILTERHDRREADSCONFIG:-}" ]]; then
    return 0
  fi

  local hdr_config_source
  hdr_config_source="$(resolve_path "${FILTERHDRREADSCONFIG}")"
  [[ -f "${hdr_config_source}" ]] || die "HDR barcode config '${hdr_config_source}' does not exist."

  HDR_FILTER_CONFIG_BASENAME="$(basename "${hdr_config_source}")"
  HDR_FILTER_CONFIG_LOCAL="${RUN_DIR}/${HDR_FILTER_CONFIG_BASENAME}"
  cp "${hdr_config_source}" "${HDR_FILTER_CONFIG_LOCAL}"
}

configure_barcode_extraction() {
  local parse_regex
  parse_regex='(?P<umi_1>[ACGT]{10})(?P<cell_3>[ACGT]{8})(?P<R2_linker>GTGGCCGATGTTTCGCATCGGCGTACGACT)(?P<cell_2>[ACGT]{8})(?P<R3_linker>ATCCACGTGCTTGAGACTGTGG)(?P<cell_1>[ACGT]{8})(?P<polyA>[ACGT]{4})(?P<polyA_2>[ACGT]{4})'

  if [[ "${LIBTYPE}" == "PARSE" ]]; then
    EXTRACT_METHOD="regex"
    BC_PATTERN="${parse_regex}"
  else
    require_var "BCPATTERN"
    EXTRACT_METHOD="string"
    BC_PATTERN="${BCPATTERN}"
  fi
}

run_whitelist_and_extract() {
  # Build the whitelist first, then extract/error-correct cell barcodes and UMI annotations into a synchronized R1/R2 output pair.
  local whitelist_args=(
    --knee-method=density
    --subset-reads=999999999
    --extract-method="${EXTRACT_METHOD}"
    --bc-pattern="${BC_PATTERN}"
    --stdin="${CLEAN_R1}"
    --plot-prefix="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}_whitelist"
    -S "${FASTQ_CORRECTED_DIR}/${SAMPLENAME}_whitelist.txt"
  )

  if [[ -n "${SETCELLNUMBER:-}" ]]; then
    whitelist_args+=(--set-cell-number="${SETCELLNUMBER}")
  else
    require_var "EXPECTEDCELLS"
    whitelist_args+=(--expect-cells="${EXPECTEDCELLS}")
  fi

  umi_tools whitelist "${whitelist_args[@]}" > "${FASTQ_CORRECTED_DIR}/${SAMPLENAME}_whitelist.log" 2>&1

  umi_tools extract \
    --extract-method="${EXTRACT_METHOD}" \
    --bc-pattern="${BC_PATTERN}" \
    --filter-cell-barcode \
    --error-correct-cell \
    --stdin="${CLEAN_R1}" \
    --stdout="${EXTRACTED_R1}" \
    --read2-in="${CLEAN_R2}" \
    --read2-out="${EXTRACTED_R2}" \
    --whitelist="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}_whitelist.txt" \
    > "${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.log" 2>&1
}

annotate_parse_barcodes_if_needed() {
  if [[ "${LIBTYPE}" != "PARSE" ]]; then
    return 0
  fi

  mv "${EXTRACTED_R2}" "${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.raw.R2.fastq.gz"
  python "${PARSE_DIR}/parse_lib_fastq_bc-to-bcID.py" \
    "${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.raw.R2.fastq.gz" \
    "${EXTRACTED_R2}" \
    "${PARSEKIT:-}"
}

split_hdr_reads_if_needed() {
  # Separate HDR-barcode reads from the CRISPResso input so downstream repair calling/crispresso alignment is not confounded by HDR barcode sequence.
  CRISPRESSO_INPUT_R2="${EXTRACTED_R2}"

  if [[ -z "${FILTERHDRREADSCONFIG:-}" ]]; then
    return 0
  fi

  CRISPRESSO_INPUT_R2="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.noBC.R2.fastq.gz"
  HDRBC_R2="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.HDRBC.R2.fastq.gz"
  HDR_BC_STATS="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.BCstats.txt"

  if [[ "${HDRBCTYPE:-}" == "cINS" ]]; then
    python "${HDR_DIR}/extractHDRreads_insBC.py" \
      -c "${HDR_FILTER_CONFIG_LOCAL}" \
      -i "${EXTRACTED_R2}" \
      -o "${CRISPRESSO_INPUT_R2}" \
      -b "${HDR_BC_STATS}" \
      -f "${HDRBC_R2}"
  elif [[ "${HDRBCTYPE:-}" == "rINS" ]]; then
    mv "${EXTRACTED_R2}" "${CRISPRESSO_INPUT_R2}"
    log "HDR barcode type is random insertion; skipping HDR barcode filtering."
  else
    python "${HDR_DIR}/extractHDRreads.py" \
      -c "${HDR_FILTER_CONFIG_LOCAL}" \
      -i "${EXTRACTED_R2}" \
      -o "${CRISPRESSO_INPUT_R2}" \
      -b "${HDR_BC_STATS}" \
      -f "${HDRBC_R2}"
  fi
}

get_crispresso_alignment_score() {
  case "${TARGET:-}" in
    GAPDH_sg13_151nt|Gpi1_200)
      printf '20\n'
      ;;
    Son_sg4|B2m_sg1|GAPDH_sg13_99nt|GAPDH_sg13_SNP|Gpi1|Pgk1)
      printf '33\n'
      ;;
    *)
      printf '41\n'
      ;;
  esac
}

run_crispresso() {
  # CRISPResso parameters vary slightly by target; centralize the branching here so the rest of the pipeline can stay linear.
  require_var "AMPLICONSEQ"
  require_var "GUIDE"
  require_var "CRISPRESSOWINDOW"

  local output_dir="${RUN_DIR}/${SAMPLENAME}_crispresso_out"
  local min_score
  min_score="$(get_crispresso_alignment_score)"

  local crispresso_args=(
    -r1 "${CRISPRESSO_INPUT_R2}"
    --amplicon_seq "${AMPLICONSEQ}"
    -o "${output_dir}"
    -g "${GUIDE}"
    -w "${CRISPRESSOWINDOW}"
    --exclude_bp_from_left 15
    --exclude_bp_from_right 0
    --plot_window_size 29
    --fastq_output
  )

  if [[ -n "${HDRSEQ:-}" ]]; then
    crispresso_args+=(-e "${HDRSEQ}")
  fi

  if [[ -n "${TARGET:-}" ]]; then
    crispresso_args+=(
      -amas "${min_score}"
      --needleman_wunsch_gap_open -40
      --needleman_wunsch_gap_incentive 2
    )
  else
    crispresso_args+=(--default_min_aln_score "${min_score}")
  fi

  CRISPResso "${crispresso_args[@]}"
}

enter_crispresso_output_dir() {
  local crispresso_dir_name
  crispresso_dir_name="$(basename "${CRISPRESSO_INPUT_R2}" .fastq.gz)"
  CRISPRESSO_RESULT_DIR="${RUN_DIR}/${SAMPLENAME}_crispresso_out/CRISPResso_on_${crispresso_dir_name}"
  cd "${CRISPRESSO_RESULT_DIR}"
}

filter_crispresso_output_if_needed() {
  if [[ "${CRISPRESSOFILTER:-}" != "GAPDH" ]]; then
    return 0
  fi

  mv CRISPResso_output.fastq.gz CRISPResso_output_unfiltered.fastq.gz
  python "${SCRIPT_DIR}/filter_crispresso_with_patterns.py" \
    "${SCRIPT_DIR}/GAPDH_pseudogenes_patterns.txt" \
    CRISPResso_output_unfiltered.fastq.gz
}

run_hdr_postprocessing() {
  # Convert HDR FASTQ output into tables, remap candidate off-target and translocation reads, then build final editing-outcome summaries.
  if [[ -z "${FILTERHDRREADSCONFIG:-}" ]]; then
    return 0
  fi

  require_var "HDRBCLENGTH"
  require_var "DROTARGET"
  require_var "DROREADLENGTH"
  require_var "CHRTARGET"
  require_var "STARTTARGET"
  require_var "ENDTARGET"
  require_var "SPECIES"

  python "${HDR_DIR}/crispresso_fastq_to_table.py" \
    "${HDRBC_R2}" \
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
    require_var "CBCPATH"
    python "${PARSE_DIR}/Parse_translocation_extractor.py" "${CBCPATH}" translocation_coordinates_filtered.csv
  else
    python "${HDR_DIR}/10X_translocation_extractor.py"
  fi

  if [[ -f unaligned_translocations.csv ]] && [[ "$(wc -l < unaligned_translocations.csv)" -gt 1 ]]; then
    SCOUT_GENOME_DIR="${GENOME_DIR}" bash "${HDR_DIR}/unaligned_translocation_bwa_mem.sh" unaligned_translocations.csv "${LIBTYPE}"
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
    require_var "CBCPATH"
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

  python "${SCRIPT_DIR}/integrate_translocations_and_filter_outcomes.py" \
    --editing-outcomes EditingOutcomesByCellBarcode.xlsx \
    --output-dir "${PWD}" \
    --target-chromosome "${CHRTARGET:-}" \
    --translocations "${translocation_summary}"
}

main() {
  local config_path="${1:-}"
  [[ -n "${config_path}" ]] || die "Usage: $0 <config-file>"
  [[ -f "${config_path}" ]] || die "Config file '${config_path}' does not exist."

  CONFIG_PATH="$(cd "$(dirname "${config_path}")" && pwd)/$(basename "${config_path}")"
  CONFIG_DIR="$(dirname "${CONFIG_PATH}")"

  source "${CONFIG_PATH}"

  require_var "SAMPLENAME"
  require_var "R1"
  require_var "R2"
  require_var "LIBTYPE"

  R1="$(resolve_path "${R1}")"
  R2="$(resolve_path "${R2}")"
  CBCPATH="$(resolve_path "${CBCPATH:-}")"
  GENOME_DIR="$(resolve_path "${GENOME_DIR:-${GENOME_DIR_DEFAULT}}")"

  [[ -f "${R1}" ]] || die "R1 FASTQ '${R1}' does not exist."
  [[ -f "${R2}" ]] || die "R2 FASTQ '${R2}' does not exist."

  RUN_DIR="${ROOT_DIR}/${SAMPLENAME}"
  if [[ -d "${RUN_DIR}" ]] && [[ "${FORCE:-FALSE}" != "TRUE" ]]; then
    log "Output directory '${RUN_DIR}' already exists. Set FORCE=TRUE to reuse it."
    exit 0
  fi

  mkdir -p "${RUN_DIR}/fastq" "${RUN_DIR}/fastq_corrected_bc"

  FASTQ_DIR="${RUN_DIR}/fastq"
  FASTQ_CORRECTED_DIR="${RUN_DIR}/fastq_corrected_bc"
  RAW_R1="${FASTQ_DIR}/${SAMPLENAME}.R1.fastq.gz"
  RAW_R2="${FASTQ_DIR}/${SAMPLENAME}.R2.fastq.gz"
  CLEAN_R1="${FASTQ_DIR}/${SAMPLENAME}.clean.R1.fastq.gz"
  CLEAN_R2="${FASTQ_DIR}/${SAMPLENAME}.clean.R2.fastq.gz"
  EXTRACTED_R1="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.R1.fastq.gz"
  EXTRACTED_R2="${FASTQ_CORRECTED_DIR}/${SAMPLENAME}.extracted.R2.fastq.gz"

  setup_logging

  ln -sfn "${R1}" "${RAW_R1}"
  ln -sfn "${R2}" "${RAW_R2}"

  prepare_hdr_filter_config

  cd "${RUN_DIR}"

  load_conda
  activate_conda_if_available "scout"
  run_seq_filter
  configure_barcode_extraction
  run_whitelist_and_extract
  annotate_parse_barcodes_if_needed
  split_hdr_reads_if_needed

  load_conda
  activate_conda_if_available "crispresso2_env"
  run_crispresso
  conda deactivate || true

  enter_crispresso_output_dir
  filter_crispresso_output_if_needed
  run_hdr_postprocessing

  log "Pipeline finished for sample '${SAMPLENAME}'."
}

main "$@"
