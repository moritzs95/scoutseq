#!/usr/bin/env bash

# Remap candidate off-target sequences back to the reference genome and
# summarize surviving coordinates per cell barcode.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/../required-core/genome_utils.sh"

if [[ "$#" -ne 6 ]]; then
  echo "Usage: $0 <chromosome> <start_position> <end_position> <input_csv> <species> <libtype>" >&2
  exit 1
fi

chromosome="$1"
start_position="$2"
end_position="$3"
input_csv="$4"
species="$5"
libtype="$6"

output_fasta="potential_offtargets.fasta.gz"
output_bam="offtargets.sorted.bam"
filtered_bam="filtered_offtargets.bam"
filtered_csv="offtargets_coordinates_filtered.csv"
genome_path="$(resolve_genome_file "${species}")"

# Convert candidate sequences to a FASTA that BWA can align.
awk -F',' 'NR>1 {print ">"$1"_"$2"_"$5"_"$6"\n"$7}' "${input_csv}" | gzip > "${output_fasta}"

bwa mem -t 24 -M -R '@RG\tID:Id\tPL:illumina\tLB:lb1\tSM:Id' "${genome_path}" "${output_fasta}" \
  | samtools view -@ 10 -bS \
  | samtools sort -@ 10 -o "${output_bam}"

samtools index "${output_bam}"

samtools view -h "${output_bam}" \
  | awk -v chr="${chromosome}" -v start="${start_position}" -v end="${end_position}" \
    '($3 == chr && $4 >= start && $4 <= end) {next} {print $0}' \
  | samtools view -bS - > "${filtered_bam}"

samtools index "${filtered_bam}"

if [[ "${libtype}" == "PARSE" ]]; then
  samtools view "${filtered_bam}" \
    | awk 'BEGIN {OFS=","; print "bc,umi,mod_offtarget,read_count,chromosome,start,end,location"} \
    {split($1, name_parts, "_"); bc=name_parts[1] "_" name_parts[2] "_" name_parts[3]; umi=name_parts[4]; mod_offtarget=name_parts[5]; read_count=name_parts[6]; chrom=$3; start=$4; end=$4+length($10)-1; location=chrom":"start"-"end; print bc,umi,mod_offtarget,read_count,chrom,start,end,location}' \
    > "${filtered_csv}"
else
  samtools view "${filtered_bam}" \
    | awk 'BEGIN {OFS=","; print "bc,umi,mod_offtarget,read_count,chromosome,start,end,location"} \
    {split($1, name_parts, "_"); bc=name_parts[1]; umi=name_parts[2]; mod_offtarget=name_parts[3]; read_count=name_parts[4]; chrom=$3; start=$4; end=$4+length($10)-1; location=chrom":"start"-"end; print bc,umi,mod_offtarget,read_count,chrom,start,end,location}' \
    > "${filtered_csv}"
fi
