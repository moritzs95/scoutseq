#!/usr/bin/python
"""Remove CRISPResso output reads that match known unwanted sequence patterns."""

import sys
import gzip
import pyfastx
import os

# Get fname from parameter
patternFile = sys.argv[1]
# Input fastq.gz file path
input_fastq_gz = sys.argv[2]

# Load patterns
patterns = set(x.strip() for x in open(patternFile))

# Function to check if a read contains any of the filtering anchors
def contains_anchor(sequence, anchors):
    for anchor in anchors:
        if anchor in sequence:
            return True
    return False



# Path to the index file
index_file_input = input_fastq_gz + '.fxi'

# Check if the index file exists and delete it if it does
if os.path.exists(index_file_input):
    os.remove(index_file_input)
    print(f"Deleted existing index file: {index_file_input}")

index_file_unfiltered = index_file_input.replace("output", "output_unfiltered")
# Check if the index file exists and delete it if it does
if os.path.exists(index_file_unfiltered):
    os.remove(index_file_unfiltered)
    print(f"Deleted existing index file: {index_file_unfiltered}")

# Read in fastq.gz file
fq_indel = pyfastx.Fastq(input_fastq_gz, build_index=True)

# Output filtered fastq.gz file path
output_fastq_gz = 'CRISPResso_output.fastq.gz'

# Open the output file for writing filtered reads
with gzip.open(output_fastq_gz, 'wt') as output_file:
    for r in fq_indel:
        lines = r.raw.split('\n')
        read_seq = lines[1]  # The sequence line
        if not contains_anchor(read_seq, patterns):
            # Write the FASTQ entry to the output file if it does not contain the anchor
            output_file.write(f"{lines[0]}\n{lines[1]}\n{lines[2]}\n{lines[3]}\n")

print(f"Filtered FASTQ saved to {output_fastq_gz}")
