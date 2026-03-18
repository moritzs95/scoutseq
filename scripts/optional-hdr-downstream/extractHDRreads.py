#!/usr/bin/env python

"""
Split substitution-style HDR barcode reads away from the CRISPResso input.
"""

import sys
import os
import gzip
import configparser
import re
from optparse import OptionParser
from collections import defaultdict

#parse 8-line concatenatedfastq
def process_fastq(lines=None):
    ks = ['name', 'seq', 'opt', 'q']
    return {k: v for k, v in zip(ks, lines)}

#calculate hamming distance between two strings
def hdist(a, b):
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(a, b))))


#anchors
#B2m    AGATAGCTGAGCAATAAATCTTCAATAAGTATTTTGATCAGAATAATAAATATAATTTTAAGAACAATAGTTGATCATATGCCAAACCCTCTGTACTTNNNNNNNNNNNNATGCAGTTACTCATCTTTGGTCTATCACAACATAAGTGACATACTTTCCTTTTGGTAAAGCAAAGAGGCCTAATTGAAGTCTGTCACTGTG
#                                                                           <------anchor------>                                <------anchor------>
# 32 bases between anchors
# 10 bases left and right of barcode
# anchor_left = AGTTGATCATATGCCAAACC (20 bases)
# anchor_right = TCATCTTTGGTCTATCACAA (20 bases)
# barcode WT = CTCATTACTTGG (12 bases)
#
#Tubb5           AAAGCAAAGTAGTGGATCAGGGATGCTAGGTAGAGACCACCAGGAAGAGAAGGGGGTGGGGTTTTCCAGTCAGGGCCATTTAGAATCCANNNNNNNNNNNNGTCAGCAGGGCTTTGTTTTGTTTTTCTCCTGCCTCATCTCTCAGCCTCAGGAGAGGTATTAACAGTATTATCTCCATTTATATCCTCCCAGCTGTCCTG
#                                                                           <------anchor------>                                <-----anchor----->
# anchor_left = GGTTTTCCAGTCAGGGCCAT
# anchor_right = CTTTGTTTTGTTTTTCTC
# barcode WT = CCTATGCTTTCA

# Function to allow up to 1 mismatch in anchor matching
def match_with_mismatch(anchor, seq, max_mismatches=1):
    # Loop over all possible substrings of the correct length in the sequence
    for i in range(len(seq) - len(anchor) + 1):
        subseq = seq[i:i + len(anchor)]
        if hdist(anchor, subseq) <= max_mismatches:
            return True, i  # Return True and the starting index of the match
    return False, -1


def checkForAnchors(seq, left_anchor, right_anchor):
    left_match, left_index = match_with_mismatch(left_anchor, seq)
    right_match, right_index = match_with_mismatch(right_anchor, seq)

    if left_match and right_match and right_index > left_index:
        pattern_seq = seq[left_index + len(left_anchor):right_index]
        return {
            'hit': True,
            'start': left_index,
            'end': right_index + len(right_anchor),
            'seq': pattern_seq,
            'len': len(pattern_seq)
        }
    else:
        return {'hit': False}


parser = OptionParser()
parser.add_option("-c", "--config", dest="config",
                  help="Config file with anchor and barcode information, see help for a template", type="str")
parser.add_option("-i", "--input", dest="input",
                  help="Input fastq file", type="str")
parser.add_option("-o", "--outfile", dest="output",
                  help="Output fastq file name", type="str")
parser.add_option("-b", "--barcodes", dest='bc',
                  help="Output barcode tab delimited file name", type="str")
parser.add_option("-f", "--barcodesfastq", dest='bcfastq',
                  help="Output barcode fastq file name", type="str")
parser.add_option("-v", "--verbose", dest="verbose",
                  help="Turn debug output on", default=False, action="store_true")
(options, args) = parser.parse_args()

#set some options manually for testing
#options.input = "/mnt/storage/markus/moritz/HDR_BC_substitutions/merged/B2m_sg1_sub_3T3.extendedFrags.fastq"
#options.output = "/mnt/storage/markus/moritz/HDR_BC_substitutions/results/B2m_sg1_sub_3T3_no_BC.fastq.gz"
#options.bc = "/mnt/storage/markus/moritz/HDR_BC_substitutions/results/B2m_sg1_sub_3T3_BC.txt"

#left_anchor = AGTTGATCATATGCCAAACC
#right_anchor = TCATCTTTGGTCTATCAC
#bp_between_anchors = 32
#bc_start = 10
#bc_length = 12
#bc_ref = CTCATTACTTGG


if options.verbose:
    print("Checking arguments")
    print(options)
    #config.get('CONFIG', 'fastq')

#Check input file path
try:
    fastqInPath = options.input
except IndexError as ie:
    raise SystemError("Error: Specify file name\n")

if not os.path.exists(fastqInPath):
    raise SystemError("Error: Input file does not exist\n")

#Check output file paths
try:
    fastqOutPath = options.output
except IndexError as ie:
    raise SystemError("Error: Specify output fastq output file name\n")
try:
    bcOutPath = options.bc
except IndexError as ie:
    raise SystemError("Error: Specify barcode output file name\n")
try:
    bcFastqOutPath = options.bcfastq
except IndexError as ie:
    raise SystemError("Error: Specify barcode fastq output file name\n")

#if not os.path.exists(os.path.dirname(fastqOutPath)):
#    os.mkdir(os.path.dirname(fastqOutPath))
#    #raise SystemError(f"Error: Output file folder {os.path.dirname(fastqOutPath)} is not accessible")

#Check config file path
try:
    conf = options.config
except IndexError as ie:
    raise SystemError("Error: Specify config file name\n")

#if not os.path.exists(conf):
#    raise SystemError("Error: Config file does not exist\n")

if options.verbose:
    print("Reading config file")
#get barcode/anchor info
config = configparser.ConfigParser()
config.read(conf)
#config.read("/mnt/storage/markus/moritz/HDR_BC_substitutions/B2m.config")

if options.verbose:
    print(f"{config.sections()}")
    print(f"{config.get('CONFIG', 'bc_ref')}")


#get config data
bp_between_anchors = int(config.get('CONFIG', 'bp_between_anchors'))
bc_start = int(config.get('CONFIG', 'bc_start'))
bc_length = int(config.get('CONFIG', 'bc_length'))
bc_ref = config.get('CONFIG', 'bc_ref')
left_anchor = config.get('CONFIG', 'left_anchor')
right_anchor = config.get('CONFIG', 'right_anchor')

fastqOut = gzip.open(fastqOutPath, 'wt')
bcFastqOut = gzip.open(bcFastqOutPath, 'wt')
bcOut = open(bcOutPath, 'w')
#write column header for anchored reads output
bcOut.write(f"Outcome\tReadName\tBarcodeSeq\tAnchor1Start\tAnchor2End\tBasesBetweenAnchors\tHdistBarcodeToWildtype\n")

n = 4
#with open(fastqInPath, 'r') as fh:
with gzip.open(fastqInPath, 'rt') as fh:
    lines = []
    for line in fh:
        lines.append(line.rstrip())
        if len(lines) == n:
            rec = process_fastq(lines)
            # print(f"{rec}")
            bcCheck = checkForAnchors(rec['seq'], left_anchor, right_anchor)
            if bcCheck['hit']:
                bc_region = bcCheck['seq'][bc_start:bc_start+bc_length]
                if bcCheck['len'] == bp_between_anchors:
                    ref_match = "WT"
                    ref_hdist = 0
                    if bc_region != bc_ref:
                        ref_match = "BC"
                        ref_hdist = hdist(bc_region, bc_ref)
                        if ref_hdist == 0:
                           ref_match = "WT"
                else:
                    ref_match = "INDEL"
                    ref_hdist = -1
                bcOut.write(f"{ref_match}\t{rec['name']}\t{bc_region}\t{bcCheck['start']}\t{bcCheck['end']}\t{bcCheck['len']}\t{ref_hdist}\n")
                if ref_match == "BC":
                    tmp_name = rec['name'].split(' ')
                    new_name = tmp_name[0]+'_'+bc_region+'_'+str(ref_hdist)+' '+tmp_name[1]
                    bcFastqOut.write(f"{new_name}\n{rec['seq']}\n{rec['opt']}\n{rec['q']}\n")
                else:
                    fastqOut.write(f"{rec['name']}\n{rec['seq']}\n{rec['opt']}\n{rec['q']}\n")
            else:
                fastqOut.write(f"{rec['name']}\n{rec['seq']}\n{rec['opt']}\n{rec['q']}\n")
            lines = []

fastqOut.close()
bcFastqOut.close()
bcOut.close()
