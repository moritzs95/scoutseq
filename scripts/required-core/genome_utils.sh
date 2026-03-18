#!/usr/bin/env bash

set -euo pipefail

resolve_genome_file() {
  local species="$1"
  local genome_dir="${SCOUT_GENOME_DIR:-$(pwd)/genome}"
  local -a candidates=()
  local candidate

  shopt -s nullglob

  case "${species}" in
    human)
      candidates=(
        "${genome_dir}"/hg38*.fa
        "${genome_dir}"/hg38*.fasta
        "${genome_dir}"/human*.fa
        "${genome_dir}"/human*.fasta
        "${genome_dir}"/genome.fa
        "${genome_dir}"/genome.fasta
      )
      ;;
    mouse)
      candidates=(
        "${genome_dir}"/mm10*.fa
        "${genome_dir}"/mm10*.fasta
        "${genome_dir}"/mouse*.fa
        "${genome_dir}"/mouse*.fasta
        "${genome_dir}"/genome.fa
        "${genome_dir}"/genome.fasta
      )
      ;;
    aav)
      candidates=(
        "${genome_dir}"/AAV*.fa
        "${genome_dir}"/AAV*.fasta
        "${genome_dir}"/aav*.fa
        "${genome_dir}"/aav*.fasta
      )
      ;;
    *)
      printf 'Unsupported species key: %s\n' "${species}" >&2
      return 1
      ;;
  esac

  for candidate in "${candidates[@]}"; do
    if [[ -f "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done

  printf 'No genome FASTA found for %s in %s\n' "${species}" "${genome_dir}" >&2
  return 1
}
