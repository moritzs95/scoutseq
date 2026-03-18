#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/../required-core/genome_utils.sh"

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 <input_csv> <libtype>" >&2
  exit 1
fi

input_csv="$1"
libtype="$2"

output_fasta="AAV_translocations.fasta.gz"
output_bam="AAV_translocations.sorted.bam"
output_csv="AAV_translocation_coordinates.csv"
genome_path="$(resolve_genome_file "aav")"

awk -F',' 'NR>1 {print ">"$1"_"$2"_"$3"_"$4"_"$5"\n"$10}' "${input_csv}" | gzip > "${output_fasta}"

bwa mem -t 24 -M -R '@RG\tID:Id\tPL:illumina\tLB:lb1\tSM:Id' "${genome_path}" "${output_fasta}" \
  | samtools view -@ 10 -bS \
  | samtools sort -@ 10 -o "${output_bam}"

samtools index "${output_bam}"

if [[ "${libtype}" == "PARSE" ]]; then
  samtools view "${output_bam}" \
    | awk 'BEGIN {OFS=","; print "bc,umi,read_seq,mod_translocation,read_count,chromosome,start,end,location,translocated_seq"} \
    {split($1, name_parts, "_"); bc=name_parts[1] "_" name_parts[2] "_" name_parts[3]; umi=name_parts[4]; read_seq=name_parts[5]; mod_translocation=name_parts[6]; read_count=name_parts[7]; chrom=$3; start=$4; end=$4+length($10)-1; location=chrom":"start"-"end; translocated_seq=$10; print bc,umi,read_seq,mod_translocation,read_count,chrom,start,end,location,translocated_seq}' \
    > "${output_csv}"
else
  samtools view "${output_bam}" \
    | awk 'BEGIN {OFS=","; print "bc,umi,read_seq,mod_translocation,read_count,chromosome,start,end,location,translocated_seq"} \
    {split($1, name_parts, "_"); bc=name_parts[1]; umi=name_parts[2]; read_seq=name_parts[3]; mod_translocation=name_parts[4]; read_count=name_parts[5]; chrom=$3; start=$4; end=$4+length($10)-1; location=chrom":"start"-"end; translocated_seq=$10; print bc,umi,read_seq,mod_translocation,read_count,chrom,start,end,location,translocated_seq}' \
    > "${output_csv}"
fi
