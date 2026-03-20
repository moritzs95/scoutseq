"""Convert HDR-enriched CRISPResso FASTQ output into tabular candidate events."""

import pyfastx
import re
from collections import defaultdict
import pandas as pd
import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "required-core"))

from bdr_utils import normalize_bdr_barcode

filename_hdr = str(sys.argv[1])
len_hdrbc = int(sys.argv[2])
libtype = sys.argv[3]
target = sys.argv[4]
read_length = int(sys.argv[5])
parsekit = sys.argv[6] if len(sys.argv) > 6 else ""

#calculate hamming distance between two strings
def hdist(a, b):
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(a, b))))

if parsekit == "mini":
    def adjust_and_reformat_barcode(barcode):
        # Split the barcode by colons to get individual parts
        parts = barcode.split(":")
        
        try:
            # Try to convert the first part to an integer
            first_part = int(parts[0])
            # If it's greater than 12, adjust it by subtracting 12
            if first_part > 12:
                parts[0] = str(first_part - 12)
        except ValueError:
            # If conversion fails, keep the first part as is
            pass
        
        # Zero-pad each part to ensure it's two digits, except for parts containing '-'
        parts = [part.zfill(2) if '-' not in part else part for part in parts]
        
        # Join the parts with underscores instead of colons
        formatted_barcode = "_".join(parts)
        
        return str(formatted_barcode)
else:
    def adjust_and_reformat_barcode(barcode):
        # Split the barcode by colons to get individual parts
        parts = barcode.split(":")
        
        try:
            # Try to convert the first part to an integer
            first_part = int(parts[0])
            # If it's greater than 96, adjust it by subtracting 96
            if first_part > 96:
                parts[0] = str(first_part - 96)
        except ValueError:
            # If conversion fails, keep the first part as is
            pass
        
        # Zero-pad each part to ensure it's two digits, except for parts containing '-'
        parts = [part.zfill(2) if '-' not in part else part for part in parts]
        
        # Join the parts with underscores instead of colons
        formatted_barcode = "_".join(parts)
        
        return str(formatted_barcode)
    
#read_seq is using the aligned sequence (not the raw sequence)

#define scOUT target (GAPDH_sg13_SNP, Son_sg4, B2m_sg1, GAPDH_sg13)
if target not in ["GAPDH_sg13_SNP", "Son_sg4", "B2m_sg1", "GAPDH_sg13", "Pgk1_sg1", "Gpi1_sg90", "Gpi1_sg90_200"]:
    sys.exit("Unrecognized scOUT target. Aborting.")
elif target == "Son_sg4":
        # Define the anchors - Son_sg4
    left_anchor = 'CACCCAGC' #PAM proximal
    right_anchor = 'CGCTATTTTG' #PAM distal
    HDR_anchor = 'GGCATGC' #PAM proximal including PAM mutation if applicable
    #Son_sg4 150 bp
    GSPrimer = "GATAGACGTAAATAAAAATGCTGTAAC"
    GSPrimer_length = len(GSPrimer)
    cdna_length = 203
    reference_length = 200
    WT_amplicon = 'GATAGACGTAAATAAAAATGCTGTAACCGACTTATCTAATAAAAATTGGCACCCAGCCGCTATTTTGTTGACTGAGGAAGTTTATGTTAATTTTTTAGGGTCTGATAGAATATTCATGTGTATTACAGTGGTATTCATATGCTATGTCTCT'
    WT_amplicon_after_cutsite = 'TTGGCACCCAGCCGCTATTTTGTTGACTGAGGAAGTTTATGTTAATTTTTTAGGGTCTGATAGAATATTCATGTGTATTACAGTGGTATTCATATGCTATGTCTCTAAACTTTATTTTCAAAAGCTTAAGGCCCAAATACAAACTTCTCTGGAAT'
    hdrbc_len = 15
    thresholds = {
    'I90plus': 30,
    'D90plus': 30,
    'S1plus': 50
    }
elif target == "B2m_sg1":
    # Define the anchors - B2m_sg1 (other direction than son)
    left_anchor = 'ACTTGGAT' #PAM proximal
    right_anchor = 'ACTTCTCATT' #PAM distal
    HDR_anchor = 'TCAATGCAGT' #PAM proximal including PAM mutation if applicable
    #B2m_sg1 150 bp
    GSPrimer = "GTATTTTGATCAGAATAATAAATATAATTTTAAGAA"
    GSPrimer_length = len(GSPrimer)
    cdna_length = 249
    reference_length = 200
    WT_amplicon = 'GTATTTTGATCAGAATAATAAATATAATTTTAAGAACAATAGTTGATCATATGCCAAACCCTCTGTACTTCTCATTACTTGGATGCAGTTACTCATCTTTGGTCTATCACAACATAAGTGACATACTTTCCTTTTGGTAAAGCAAAGAGGC'
    WT_amplicon_after_cutsite = 'CCTCTGTACTTCTCATTACTTGGATGCAGTTACTCATCTTTGGTCTATCACAACATAAGTGACATACTTTCCTTTTGGTAAAGCAAAGAGGCCTAATTGAAGTCTGTCACTGTGCCCAATGCTTAGCAATTCTCACCCCCA'
    hdrbc_len = 15
    thresholds = {
    'I90plus': 30,
    'D90plus': 30,
    'S1plus': 50
    }
elif target == "Pgk1_sg1":
    # Define the anchors - B2m_sg1 (other direction than son)
    left_anchor = 'GGTTCCTGTG' #PAM proximal
    right_anchor = 'CTCCTAAGT' #PAM distal
    HDR_anchor = 'GGTTTTAGTG' #not relevant
    #Pgk1_sg1 91 bp
    GSPrimer = "TCCTTCCTGGGGTGGATGCT"
    GSPrimer_length = len(GSPrimer)
    cdna_length = 455
    reference_length = 140
    WT_amplicon = 'TCCTTCCTGGGGTGGATGCTCTCAGCAATGTTTAGTATTTTCTTTCCTGCCTTTGGTTCCTGTGCTCCTAAGTCAACCTAGTGTTTTCCACATCTCCATTTGGTGTTAGCGCAAGATTCAGCTAGTGGCTGAGATGTGGC'
    WT_amplicon_after_cutsite = 'TTGGTTCCTGTGCTCCTAAGTCAACCTAGTGTTTTCCACATCTCCATTTGGTGTTAGCGCAAGATTCAGCTAGTGGCTGAGATGTGGC'
    hdrbc_len = 12
    thresholds = {
    'I90plus': 30,
    'D90plus': 30,
    'S1plus': 50
    }
elif target == "Gpi1_sg90":
    # Define the anchors - B2m_sg1 (other direction than son)
    left_anchor = 'GTCCTCCG' #PAM proximal
    right_anchor = 'TGTCCCTTCT' #PAM distal
    HDR_anchor = 'GTAAGCCG' #not relevant
    #Gpi1_sg90 91 bp
    GSPrimer = "AAGCAACAGCGGGACACCAA"
    GSPrimer_length = len(GSPrimer)
    cdna_length = 1114
    reference_length = 140
    WT_amplicon = 'AAGCAACAGCGGGACACCAAACTAGAATAACTCCAGCCGCGGCCCTACTGACTGGTCCTCCGTGTCCCTTCTCACCATATGCACTGCATGGTCCTGCCCCTCCCTGCCCAGAGCGCACCACCGGTAGTTGGCCTGGACTA'
    WT_amplicon_after_cutsite = 'TGTCCCTTCTCACCATATGCACTGCATGGTCCTGCCCCTCCCTGCCCAGAGCGCACCACCGGTAGTTGGCCTGGACTA'
    hdrbc_len = 12
    thresholds = {
    'I90plus': 30,
    'D90plus': 30,
    'S1plus': 50
    }
elif target == "Gpi1_sg90_200":
    # Define the anchors 
    left_anchor = 'GTCCTCCG' #PAM proximal
    right_anchor = 'TGTCCCTTCT' #PAM distal
    HDR_anchor = 'GTAAGCCG' #relevant!
    #Gpi1_sg90 91 bp
    GSPrimer = "AAGCAACAGCGGGACACCAA"
    GSPrimer_length = len(GSPrimer)
    cdna_length = 1114
    reference_length = 200
    WT_amplicon = 'AAGCAACAGCGGGACACCAAACTAGAATAACTCCAGCCGCGGCCCTACTGACTGGTCCTCCGTGTCCCTTCTCACCATATGCACTGCATGGTCCTGCCCCTCCCTGCCCAGAGCGCACCACCGGTAGTTGGCCTGGACTACAAGGCTGTTGGGAGAAGCTGGTCTGGAACTGCCATCCACCCACTACGCACCCTCCCTGT'
    WT_amplicon_after_cutsite = 'TGTCCCTTCTCACCATATGCACTGCATGGTCCTGCCCCTCCCTGCCCAGAGCGCACCACCGGTAGTTGGCCTGGACTACAAGGCTGTTGGGAGAAGCTGGTCTGGAACTGCCATCCACCCACTACGCACCCTCCCTGT'
    hdrbc_len = 15
    thresholds = {
    'I90plus': 30,
    'D90plus': 30,
    'S1plus': 50
    }
elif target == "GAPDH_sg13":
    # Define the anchors - GAPDH_sg13 (other direction than son)
    left_anchor = 'CTCCTCAC' #PAM proximal
    right_anchor = 'AGTTGCCATG' #PAM distal
    HDR_anchor = 'CTCCCCTTGT' #PAM proximal including PAM mutation
    #GAPDH sg13 151 bp
    GSPrimer = "GTCCCTGCCACACTCAGTCCCCC"
    GSPrimer_length = len(GSPrimer)
    cdna_length = 137
    reference_length = 200
    WT_amplicon = 'GTCCCTGCCACACTCAGTCCCCCACCACACTGAATCTCCCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGTACCATCAATAAAGTACCCTGTGCTCAACCAGTTACTTGTCCTGT'
    WT_amplicon_after_cutsite = 'CCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGTACCATCAATAAAGTACCCTGTGCTCAACCAGTTACTTGTCCTGTCTTATTCTAGGGTCTGGGGCAGAGGGGAGGGAAGCTGGGCTTGTGTCAA'
    hdrbc_len = 15
    thresholds = {
    'I90plus': 30,
    'D90plus': 30,
    'S1plus': 50
    }
    if read_length == 99 or read_length == 90:
        cdna_length = 137
        reference_length = 137
        thresholds = {
        'I90plus': 30,
        'D90plus': 30,
        'S1plus': 50
        }
elif target == "GAPDH_sg13_SNP":
    #GAPDH_sg13_SNP
    left_anchor = 'CTCCTCAC'  # PAM proximal
    right_anchor = 'AGTTTCCATG'  # PAM distal
    HDR_anchor = 'CTCCCCTACT'  # PAM proximal including PAM mutation
    #GAPDH SNP 102 bp
    GSPrimer = "GTCCCTGCCACACTCAGTCCCCC"
    GSPrimer_length = len(GSPrimer)
    cdna_length = 137
    reference_length = 137
    WT_amplicon = 'GTCCCTGCCACACTCAGTCCCCCACCACACTGAATCTCCCCTCCTCACAGTTTCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGT'
    WT_amplicon_after_cutsite = 'CCCTCCTCACAGTTTCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGTACCATCAATAAAGTACCCTGTGCTCAACCAGTTACTTGTCCTGTCTTATTCTAGGGTCTGGGGCAGAGGGGAGGGAAGCTGGGCTTGTGTCAA'
    hdrbc_len = 12
    thresholds = {
    'I90plus': 30,
    'D90plus': 30,
    'S1plus': 50
    }


# Function to check if the row should be filtered out based only on the first mod entry
#filter deletions with only trailing gaps
def all_gaps_trailing(seq, mod):
    if mod.startswith('D'):
        # Find the first occurrence of a gap
        first_gap_index = seq.find('-')
        
        # If no gap is found, return True (since there are no gaps)
        if first_gap_index == -1:
            return True
        
        # Check if all characters after the first gap are gaps
        return all(char == '-' for char in seq[first_gap_index:])
    else:
        return False
    

# General function to calculate deletion thresholds
def calculate_threshold_deletions(number, read_length, cdna_length, reference_length):
    effective_read_length = min(read_length + number, cdna_length, reference_length)
    match = effective_read_length - number
    tolerance = 5
    aln_score_threshold = (match - tolerance) / reference_length * 100
    return max(aln_score_threshold, 24)

# General function to calculate insertion thresholds
def calculate_threshold_insertions(number, read_length, cdna_length, reference_length):
    effective_read_length = min(cdna_length + number, read_length)
    match = effective_read_length - number
    tolerance = 5
    aln_score_threshold = (match - tolerance) / (reference_length + number) * 100
    return max(aln_score_threshold, 30)

# Function to generate thresholds for a given target
def generate_thresholds(read_length, cdna_length, reference_length):
    thresholds = {}
    for number in range(1, 91):  # Generate thresholds for D1 to D90 and I1 to I90
        key_del = f'D{number}'
        key_ins = f'I{number}'
        thresholds[key_del] = calculate_threshold_deletions(number, read_length, cdna_length, reference_length)
        thresholds[key_ins] = calculate_threshold_insertions(number, read_length, cdna_length, reference_length)
    
    return thresholds

new_thresholds = generate_thresholds(read_length, cdna_length, reference_length)
thresholds.update(new_thresholds)

print(thresholds)

if os.path.exists(filename_hdr):
    # Path to the index file
    index_file_hdr = filename_hdr + '.fxi'

    # Check if the index file exists and delete it if it does
    if os.path.exists(index_file_hdr):
        os.remove(index_file_hdr)
        print(f"Deleted existing index file: {index_file_hdr}")

    #read in hdr barcode read data
    fq_hdr = pyfastx.Fastq(filename_hdr, build_index=True)
    #create nested dict
    bc_dict = defaultdict(lambda: defaultdict())

    if parsekit == "mini":
        def adjust_and_reformat_barcode(barcode):
            # Split the barcode by colons to get individual parts
            parts = barcode.split(":")
            
            try:
                # Try to convert the first part to an integer
                first_part = int(parts[0])
                # If it's greater than 12, adjust it by subtracting 12
                if first_part > 12:
                    parts[0] = str(first_part - 12)
            except ValueError:
                # If conversion fails, keep the first part as is
                pass
            
            # Zero-pad each part to ensure it's two digits, except for parts containing '-'
            parts = [part.zfill(2) if '-' not in part else part for part in parts]
            
            # Join the parts with underscores instead of colons
            formatted_barcode = "_".join(parts)
            
            return str(formatted_barcode)
    else:
        def adjust_and_reformat_barcode(barcode):
            # Split the barcode by colons to get individual parts
            parts = barcode.split(":")
            
            try:
                # Try to convert the first part to an integer
                first_part = int(parts[0])
                # If it's greater than 96, adjust it by subtracting 96
                if first_part > 96:
                    parts[0] = str(first_part - 96)
            except ValueError:
                # If conversion fails, keep the first part as is
                pass
            
            # Zero-pad each part to ensure it's two digits, except for parts containing '-'
            parts = [part.zfill(2) if '-' not in part else part for part in parts]
            
            # Join the parts with underscores instead of colons
            formatted_barcode = "_".join(parts)
            
            return str(formatted_barcode)



    for r in fq_hdr:
        lines = r.raw.split('\n')
        parts = lines[0].split()[0].split('_')
        if libtype == "10X":
            name,bc,umi,hdrbc,hd = parts
        elif libtype == "BDR":
            name,bc,umi,hdrbc,hd = parts
            raw_bc = bc
            bc = normalize_bdr_barcode(bc)
        else:
            name,bc,bcSEQ,stype,umi,hdrbc,hd = parts
            bc = adjust_and_reformat_barcode(bc)
            umi = str(umi)+"-"+str(stype)
        hdr_bc = str(hdrbc)
        raw_read_seq = lines[1]
        if not umi in bc_dict[bc]:
            bc_dict[bc][umi] = defaultdict(int)
        if libtype == "BDR":
            bc_dict[bc][umi]["raw_bc"] = raw_bc
        if len(hdrbc) == len_hdrbc:
            if not 'hdr_bc' in bc_dict[bc][umi]:
                bc_dict[bc][umi]['hdr_bc'] = defaultdict(int)
            if not hdr_bc in bc_dict[bc][umi]['hdr_bc']:
                bc_dict[bc][umi]['hdr_bc'][hdr_bc] = {'count': 1, raw_read_seq: 1, 'hd': int(hd)}
            else:
                bc_dict[bc][umi]['hdr_bc'][hdr_bc]['count'] += 1
                if raw_read_seq in bc_dict[bc][umi]['hdr_bc'][hdr_bc].keys():
                    bc_dict[bc][umi]['hdr_bc'][hdr_bc][raw_read_seq] += 1
                else:
                    bc_dict[bc][umi]['hdr_bc'][hdr_bc][raw_read_seq] = 1



    for bc in bc_dict:
        for umi in bc_dict[bc]:
            if 'hdr_bc' in bc_dict[bc][umi]:
                max_hdr_bc = ''
                max_hdr_bc_count = 0
                max_readseq = ""
                max_readseq_count = 0
                max_hd = 0

                # Identify the hdr_bc with the highest read count
                for hdr_bc, inner_dict in bc_dict[bc][umi]['hdr_bc'].items():
                    current_hdr_bc_count = inner_dict['count']
                    current_max_readseq_count = 0
                    current_max_readseq = ""
                    
                    # Check each raw_read_seq and its count to find the most common one
                    for key, value in inner_dict.items():
                        if key not in ("count", "hd") and value > current_max_readseq_count:
                            current_max_readseq_count = value
                            current_max_readseq = key
                    
                    # Check if this hdr_bc has a higher read count than the current max
                    if current_max_readseq_count > max_readseq_count:
                        max_hdr_bc = hdr_bc
                        max_hdr_bc_count = current_hdr_bc_count
                        max_readseq = current_max_readseq
                        max_readseq_count = current_max_readseq_count
                        max_hd = inner_dict['hd']

                # Overwrite bc_dict[bc][umi]['hdr_bc'] with only the maximum entry
                new_hdr_bc = max_hdr_bc + "_" + str(max_hdr_bc_count) + "_" + max_readseq + "_" + str(max_readseq_count) + "_" + str(max_hd)
                
                # Instead of keeping all hdr_bc entries, retain only the one with max count
                bc_dict[bc][umi]['hdr_bc'] = new_hdr_bc

    
    bcs = []
    frames = []
    for bc, d in bc_dict.items():
        bcs.append(bc)
        frames.append(pd.DataFrame.from_dict(d, orient='index'))


    dfmain_hdr = pd.concat(frames, keys=bcs).reset_index()

    # Rename the index columns that were introduced by concat/reset_index.
    rename_map = {}
    if "level_0" in dfmain_hdr.columns:
        rename_map["level_0"] = "bc"
    if "level_1" in dfmain_hdr.columns:
        rename_map["level_1"] = "umi"
    dfmain_hdr = dfmain_hdr.rename(columns=rename_map)
        
    # Save the result to CSV
    dfmain_hdr.to_csv('bc_umi_hdr_seq.csv', index=False)


    

# Path to the index file
filename_indel = 'CRISPResso_output.fastq.gz'
index_file_indel = filename_indel + '.fxi'

# Check if the index file exists and delete it if it does
if os.path.exists(index_file_indel):
    os.remove(index_file_indel)
    print(f"Deleted existing index file: {index_file_indel}")



#read in indel read data
fq_indel = pyfastx.Fastq(filename_indel, build_index=True)



#for deletions
def extract_bases_after_gap(aln_seq, num_bases=15):
    # Find the position of the first gap
    first_gap_pos = aln_seq.find('-')

    # Remove the gaps and extract the first 15 bases after the first gap
    clean_seq = aln_seq[first_gap_pos:].replace('-', '')  # Remove gaps after the first gap
    return clean_seq[:num_bases]  # Extract first 'num_bases' bases

#for deletion + substitution
def extract_last_20bp_and_first_10bp(aln_seq):
    # Remove all gaps, including trailing gaps
    clean_seq = aln_seq.replace('-', '')

    # Extract the last 20 bases and then the first 10 bases of those
    return clean_seq[-20:][:10]

potential_translocations = defaultdict(lambda: defaultdict(lambda: {
    'count': 0,
    'read_seqs': defaultdict(int),
    'mods': defaultdict(int)
}))

bc_dict = defaultdict(lambda: defaultdict(lambda: {
    'aln_scores_sum': 0,
    'aln_scores_count': 0,
    'Reference_MODIFIED': 0,
    'Reference_UNMODIFIED': 0,
    'mods': None,
    'wt_seq': defaultdict(lambda: {'count': 0, 'read_seqs': defaultdict(int)})
}))

for r in fq_indel:
    lines = r.raw.split('\n')
    parts = lines[0].split()[0].split('_')
    
    if libtype == "10X":
        name, bc, umi = parts
    elif libtype == "BDR":
        name, bc, umi = parts
        raw_bc = bc
        bc = normalize_bdr_barcode(bc)
    else:
        name, bc, bcSEQ, stype, umi = parts
        bc = adjust_and_reformat_barcode(bc)
        umi = str(umi)+"-"+str(stype)
    
    raw_read_seq = lines[1]
    readlength = len(raw_read_seq)
    
    
    crispresso = dict(re.findall(r'(\S+)=(".*?"|\S+)', lines[2]))
    
    if 'ALN_SEQ' in crispresso:
        read_seq = str(crispresso['ALN_SEQ'])
        aln_ref_seq = str(crispresso['ALN_REF'])
        # Check for perfect alignment of primer sequence to catch off-targets
        seq_to_match = raw_read_seq[:GSPrimer_length]
        
        # Continue if the hamming distance is greater than 1 (off-target priming read)
        if hdist(GSPrimer, seq_to_match) > 1:
            continue

    
    if 'ALN_SCORES' in crispresso:
        aln_score = float(crispresso['ALN_SCORES'])
        
        # Only increment aln_scores_count if the CLASS is not Reference_UNMODIFIED with ALN_SCORES < 60
        if 'CLASS' in crispresso:
            if crispresso['CLASS'] == 'Reference_UNMODIFIED' and aln_score < 55:
                continue  # Skip further processing for unmodified reads with a low alignment score
            else:
                bc_dict[bc][umi]['aln_scores_count'] += 1  # Only count valid entries
                bc_dict[bc][umi]['aln_scores_sum'] += aln_score
        # If CLASS is missing, skip further processing
        else:
            continue

    if libtype == "BDR":
        bc_dict[bc][umi]["raw_bc"] = raw_bc
    
    if 'CLASS' in crispresso:
        if crispresso['CLASS'] == 'Reference_UNMODIFIED' and float(crispresso['ALN_SCORES']) < 55:
            continue  # Skip further processing if CLASS is Reference_UNMODIFIED and ALN_SCORES < 55
        
        # Increment or initialize the count for CLASS
        bc_dict[bc][umi][crispresso['CLASS']] += 1
        
        # If the class is Reference_UNMODIFIED, add or increment the wild-type sequence
        if crispresso['CLASS'] == 'Reference_UNMODIFIED':
            if bc_dict[bc][umi]['wt_seq'] is None:
                bc_dict[bc][umi]['wt_seq'] = defaultdict(lambda: {'count': 0, 'read_seqs': defaultdict(int)})
    
            bc_dict[bc][umi]['wt_seq'][read_seq]['count'] += 1
            bc_dict[bc][umi]['wt_seq'][read_seq]['read_seqs'][read_seq] += 1
    
    if 'MODS' in crispresso:
        # Find all modifications (D, I, S) and their respective counts
        combined_edits = dict(re.findall(r'([DIS])([0-9]+)', crispresso['MODS']))
        
        mod_d = 'D' + combined_edits.get('D', '0')
        mod_i = 'I' + combined_edits.get('I', '0')
        mod_s = 'S' + combined_edits.get('S', '0')

        # Initialize the mods string for this read
        mod_str = ''

        # Create a dictionary for all mods and their numeric values
        mod_values = {
            'D': int(mod_d[1:]),
            'I': int(mod_i[1:]),
            'S': int(mod_s[1:])
        }

        # Extract only non-zero mods for processing
        non_zero_mods = {mod: value for mod, value in mod_values.items() if value != 0}

        # Combine the mods into a single string for easier processing
        mod_str_temp = ''.join([f"{mod}{value}" for mod, value in non_zero_mods.items()])
        if len(non_zero_mods) == 1:
            # Process single modification (e.g., D5, I4, S1)
            mod_value = 0
            if mod_d != 'D0' and mod_i == 'I0':
                #deletions
                mod_value = int(mod_d[1:])
                effective_read_length = min(read_length + int(mod_d[1:]), cdna_length, reference_length)
                trimmed_length = max(effective_read_length - readlength, 0)
                # Process Deletion based on alignment score thresholds
                if all_gaps_trailing(read_seq, mod_d):
                    bc_dict[bc][umi]['aln_scores_count'] -= 1
                    bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                    bc_dict[bc][umi][crispresso['CLASS']] -= 1
                    continue
                elif mod_value < 90:
                    matches = thresholds[mod_d] / 100 * reference_length
                    adjusted_threshold = (matches - trimmed_length) / reference_length * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        continue
                    else:
                        extracted_bases = extract_bases_after_gap(read_seq)
                        if extracted_bases in WT_amplicon_after_cutsite:
                            mod_str = mod_str_temp
                        else:
                            bc_dict[bc][umi]['aln_scores_count'] -= 1
                            bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                            bc_dict[bc][umi][crispresso['CLASS']] -= 1
                            continue
                else:
                    matches = thresholds['D90plus'] / 100 * reference_length
                    adjusted_threshold = (matches - trimmed_length) / reference_length * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        continue
                    else:
                        extracted_bases = extract_bases_after_gap(read_seq)
                        if extracted_bases in WT_amplicon_after_cutsite:
                            mod_str = mod_str_temp
                        else:
                            bc_dict[bc][umi]['aln_scores_count'] -= 1
                            bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                            bc_dict[bc][umi][crispresso['CLASS']] -= 1
                            continue

            if mod_i != 'I0':
                mod_value = int(mod_i[1:])
                effective_read_length = min(cdna_length + int(mod_i[1:]), read_length)
                trimmed_length = max(effective_read_length - readlength, 0)
                # Process Insertion based on alignment score thresholds
                if mod_value < 90:
                    if aln_score < thresholds[mod_i]:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        continue
                    else:
                        mod_str = mod_str_temp
                else:
                    if aln_score < thresholds['I90plus']:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        continue
                    else:
                        mod_str = mod_str_temp

            if mod_s != 'S0':
                mod_value = int(mod_s[1:])
                # Process Substitution (no specific threshold handling)
                if aln_score < 55:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        continue
                else:
                    mod_str = mod_str_temp

        elif len(non_zero_mods) > 1:
            # Process complex modification (e.g., D5I4S1)
            # Combine the mods into a single string for easier processing
            if 'D' in non_zero_mods and 'S' in non_zero_mods and 'I' not in non_zero_mods:
                # Case: D and S are non-zero, I is zero
                if int(mod_d[1:]) < 90:
                    effective_read_length = min(read_length + int(mod_d[1:]), cdna_length, reference_length)
                    trimmed_length = max(effective_read_length - readlength, 0)
                    matches = thresholds[mod_d] / 100 * reference_length
                    adjusted_threshold = (matches - int(mod_s[1:]) - trimmed_length) / reference_length * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        potential_translocations[bc][umi]['count'] += 1
                        potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                        potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}S{mod_s[1:]}"] += 1
                        continue
                    else:
                        extracted_bases = extract_last_20bp_and_first_10bp(read_seq)
                        if extracted_bases in WT_amplicon_after_cutsite:
                            mod_str = mod_str_temp  # Keep the name the same (D1S1)
                        else:
                            bc_dict[bc][umi]['aln_scores_count'] -= 1
                            bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                            bc_dict[bc][umi][crispresso['CLASS']] -= 1
                            potential_translocations[bc][umi]['count'] += 1
                            potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                            potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}S{mod_s[1:]}"] += 1
                            continue
                else:
                    effective_read_length = min(read_length + int(mod_d[1:]), cdna_length, reference_length)
                    trimmed_length = max(effective_read_length - readlength, 0)
                    matches = thresholds['D90plus'] / 100 * reference_length
                    adjusted_threshold = (matches - int(mod_s[1:]) - trimmed_length) / reference_length * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        potential_translocations[bc][umi]['count'] += 1
                        potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                        potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}S{mod_s[1:]}"] += 1
                        continue
                    else:
                        extracted_bases = extract_last_20bp_and_first_10bp(read_seq)
                        if extracted_bases in WT_amplicon_after_cutsite:
                            mod_str = mod_str_temp  # Keep the name the same (D1S1)
                        else:
                            bc_dict[bc][umi]['aln_scores_count'] -= 1
                            bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                            bc_dict[bc][umi][crispresso['CLASS']] -= 1
                            potential_translocations[bc][umi]['count'] += 1
                            potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                            potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}S{mod_s[1:]}"] += 1
                            continue
            elif 'I' in non_zero_mods and 'S' in non_zero_mods and 'D' not in non_zero_mods:
                # Case: I and S are non-zero, D is zero
                if int(mod_i[1:]) < 90:
                    effective_read_length = min(cdna_length + int(mod_i[1:]), read_length)
                    trimmed_length = max(effective_read_length - readlength, 0)
                    matches = thresholds[mod_i] / 100 * (reference_length + int(mod_i[1:]))
                    adjusted_threshold = (matches - int(mod_s[1:]) - trimmed_length) / (reference_length + int(mod_i[1:])) * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        potential_translocations[bc][umi]['count'] += 1
                        potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                        potential_translocations[bc][umi]['mods'][f"I{mod_i[1:]}S{mod_s[1:]}"] += 1
                        continue
                    else:
                        mod_str = mod_str_temp  # Keep the name the same (I1S1)
                else:
                    effective_read_length = min(cdna_length + int(mod_i[1:]), read_length)
                    trimmed_length = max(effective_read_length - readlength, 0)
                    matches = thresholds['I90plus'] / 100 * (reference_length + int(mod_i[1:]))
                    adjusted_threshold = (matches - int(mod_s[1:]) - trimmed_length) / (reference_length + int(mod_i[1:])) * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        potential_translocations[bc][umi]['count'] += 1
                        potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                        potential_translocations[bc][umi]['mods'][f"I{mod_i[1:]}S{mod_s[1:]}"] += 1
                        continue
                    else:
                        mod_str = mod_str_temp  # Keep the name the same (I1S1)
            elif 'D' in non_zero_mods and 'I' in non_zero_mods and 'S' not in non_zero_mods:
                # Case: D and I are non-zero, S is zero
                if int(mod_i[1:]) < 90:
                    effective_read_length = min(cdna_length + int(mod_i[1:]), read_length)
                    trimmed_length = max(effective_read_length - readlength, 0)
                    matches = thresholds[mod_i] / 100 * (reference_length + int(mod_i[1:]))
                    adjusted_threshold = (matches - int(mod_d[1:]) - trimmed_length) / (reference_length + int(mod_i[1:])) * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        potential_translocations[bc][umi]['count'] += 1
                        potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                        potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}I{mod_i[1:]}"] += 1
                        continue
                    else:
                        mod_str = mod_str_temp # Keep the name the same 
                else:
                    effective_read_length = min(cdna_length + int(mod_i[1:]), read_length)
                    trimmed_length = max(effective_read_length - readlength, 0)
                    matches = thresholds['I90plus'] / 100 * (reference_length + int(mod_i[1:]))
                    adjusted_threshold = (matches - int(mod_d[1:]) - trimmed_length) / (reference_length + int(mod_i[1:])) * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        potential_translocations[bc][umi]['count'] += 1
                        potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                        potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}I{mod_i[1:]}"] += 1
                        continue
                    else:
                        mod_str = mod_str_temp # Keep the name the same (I1S1)
            elif 'D' in non_zero_mods and 'I' in non_zero_mods and 'S' in non_zero_mods:
                # Case: D, I, and S are all non-zero
                if int(mod_i[1:]) < 90:
                    effective_read_length = min(cdna_length + int(mod_i[1:]), read_length)
                    trimmed_length = max(effective_read_length - readlength, 0)
                    matches = thresholds[mod_i] / 100 * (reference_length + int(mod_i[1:]))
                    adjusted_threshold = (matches - int(mod_d[1:]) - int(mod_s[1:]) - trimmed_length) / (reference_length + int(mod_i[1:])) * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        potential_translocations[bc][umi]['count'] += 1
                        potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                        potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}I{mod_i[1:]}"] += 1
                        continue
                    else:
                        extracted_bases = extract_last_20bp_and_first_10bp(read_seq)
                        if extracted_bases in WT_amplicon_after_cutsite:
                            mod_str = mod_str_temp # Keep the name the same (D1I1S1)
                        else:
                            bc_dict[bc][umi]['aln_scores_count'] -= 1
                            bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                            bc_dict[bc][umi][crispresso['CLASS']] -= 1
                            potential_translocations[bc][umi]['count'] += 1
                            potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                            potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}I{mod_i[1:]}S{mod_s[1:]}"] += 1
                            continue
                        
                else:
                    effective_read_length = min(cdna_length + int(mod_i[1:]), read_length)
                    trimmed_length = max(effective_read_length - readlength, 0)
                    matches = thresholds['I90plus'] / 100 * (reference_length + int(mod_i[1:]))
                    adjusted_threshold = (matches - int(mod_d[1:]) - int(mod_s[1:]) - trimmed_length) / (reference_length + int(mod_i[1:])) * 100
                    if aln_score < adjusted_threshold:
                        bc_dict[bc][umi]['aln_scores_count'] -= 1
                        bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                        bc_dict[bc][umi][crispresso['CLASS']] -= 1
                        potential_translocations[bc][umi]['count'] += 1
                        potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                        potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}I{mod_i[1:]}"] += 1
                        continue
                    else:
                        extracted_bases = extract_last_20bp_and_first_10bp(read_seq)
                        if extracted_bases in WT_amplicon_after_cutsite:
                            mod_str = mod_str_temp # Keep the name the same (D1I1S1)
                        else:
                            bc_dict[bc][umi]['aln_scores_count'] -= 1
                            bc_dict[bc][umi]['aln_scores_sum'] -= aln_score
                            bc_dict[bc][umi][crispresso['CLASS']] -= 1
                            potential_translocations[bc][umi]['count'] += 1
                            potential_translocations[bc][umi]['read_seqs'][raw_read_seq] += 1
                            potential_translocations[bc][umi]['mods'][f"D{mod_d[1:]}I{mod_i[1:]}S{mod_s[1:]}"] += 1
                            continue
            else:
                print(mod_str)

        # Only process if there are any valid mods (i.e., mod_str is not empty)
        if mod_str:
            if bc_dict[bc][umi]['mods'] is None:
                bc_dict[bc][umi]['mods'] = defaultdict(lambda: {'count': 0, 'read_seqs': defaultdict(int), 'aln_ref_seq': defaultdict(int)})
            
            bc_dict[bc][umi]['mods'][mod_str]['count'] += 1
            bc_dict[bc][umi]['mods'][mod_str]['read_seqs'][read_seq] += 1
            bc_dict[bc][umi]['mods'][mod_str]['aln_ref_seq'][aln_ref_seq] += 1

# Processing the mods, wt_seq, and calculating average alignment score
for bc in bc_dict:
    for umi in bc_dict[bc]:
        # Process the 'mods'
        if bc_dict[bc][umi]['mods']:
            max_mods = None
            max_mods_count = 0
            max_readseq = None
            max_readseq_count = 0
            max_alnrefseq = None

            # Find the mod entry with the highest count
            for mods, inner_dict in bc_dict[bc][umi]['mods'].items():
                mods_count = inner_dict['count']
                readseq_count = max(inner_dict['read_seqs'].values())
                readseq = max(inner_dict['read_seqs'], key=inner_dict['read_seqs'].get)
                alnrefseq_count = max(inner_dict['aln_ref_seq'].values())
                alnrefseq = max(inner_dict['aln_ref_seq'], key=inner_dict['aln_ref_seq'].get)

                if mods_count > max_mods_count:
                    max_mods_count = mods_count
                    max_mods = mods
                    max_readseq = readseq
                    max_readseq_count = readseq_count
                    max_alnrefseq = alnrefseq

            # Set 'mods' field to the mod with the highest count
            if max_mods:
                bc_dict[bc][umi]['mods'] = f"{max_mods}_{max_mods_count}_{max_readseq}_{max_readseq_count}_{max_alnrefseq}"
                
        # Process the 'wt_seq' in the same manner
        # Process 'wt_seq' to keep only the sequence with the highest count
        if bc_dict[bc][umi]['wt_seq']:
            max_wtseq = None
            max_wtseq_count = 0
            for wt_seq, inner_dict in bc_dict[bc][umi]['wt_seq'].items():
                wt_count = inner_dict['count']
                if wt_count > max_wtseq_count:
                    max_wtseq_count = wt_count
                    max_wtseq = wt_seq
            
            # Set the 'wt_seq' field to the sequence with the highest count
            if max_wtseq:
                bc_dict[bc][umi]['wt_seq'] = f"{max_wtseq}_{max_wtseq_count}"
        

        # Calculate and store the average alignment score for the umi
        if bc_dict[bc][umi]['aln_scores_count'] > 0:
            avg_aln_score = bc_dict[bc][umi]['aln_scores_sum'] / bc_dict[bc][umi]['aln_scores_count']
            bc_dict[bc][umi]['avg_aln_score'] = avg_aln_score
            
# Convert the dictionary to a DataFrame
keys_to_check = ['aln_scores_sum', 'aln_scores_count', 'Reference_MODIFIED', 'Reference_UNMODIFIED']


rows = []
for bc, umi_data in bc_dict.items():
    for umi, values in umi_data.items():
        # Check if all specified keys are 0
        if all(values[key] == 0 for key in keys_to_check):
            continue  # Skip this row if all values are 0

        # Otherwise, append the row
        row = {
            'bc': bc,
            'umi': umi,
            'aln_scores_sum': values['aln_scores_sum'],
            'aln_scores_count': values['aln_scores_count'],
            'Reference_MODIFIED': values['Reference_MODIFIED'],
            'Reference_UNMODIFIED': values['Reference_UNMODIFIED'],
            'mods': values['mods'],
            'wt_seq': values['wt_seq'],  # Include wild-type sequence counts
            'avg_aln_score': values.get('avg_aln_score', 0),
            'raw_bc': values.get('raw_bc', ""),
        }
        rows.append(row)

# Function to find the index of the first mismatch between two sequences
def find_first_mismatch(seq1, seq2):
    min_length = min(len(seq1), len(seq2))  # Compare up to the length of the shorter sequence
    for i in range(min_length):
        if seq1[i] != seq2[i]:
            return i  # Return the index of the first mismatch
    return None  # Return None if no mismatch is found within the shorter sequence


# Process and count potential_translocations
translocation_rows = []
for bc, umi_data in potential_translocations.items():
    for umi, values in umi_data.items():
        # Initialize variables to track the highest count
        highest_count = 0
        best_read_seq = None
        best_mod = None
        best_mod_count = 0
        best_translocated_seq = None

        # Loop through all read sequences for the current bc and umi
        for read_seq, count in values['read_seqs'].items():
            if count > highest_count:
                # Update the best read_seq, mod, and counts if current count is higher
                highest_count = count
                first_mismatch_index = find_first_mismatch(read_seq, WT_amplicon)
                best_translocated_seq = read_seq[first_mismatch_index:]
                
                # Get the first mod and its count efficiently
                if values['mods']:
                    best_mod, best_mod_count = next(iter(values['mods'].items()))
                else:
                    best_mod, best_mod_count = "None", 0
                best_read_seq = read_seq

        # After the loop, add the best translocation information to the rows
        if best_read_seq:
            translocation_rows.append({
                'bc': bc,
                'umi': umi,
                'read_seq': best_read_seq,
                'count': highest_count,
                'mods': best_mod,
                'mods_count': best_mod_count,
                'translocated_seq': best_translocated_seq
            })

#differentiate translocations and off target amplifications
def classify_translocation(row):
    read_len = len(row['read_seq'])
    translocated_len = len(row['translocated_seq'])
    
    # Calculate the fraction
    fraction = translocated_len / read_len
    
    # Return "translocation" if the fraction is < 0.7, otherwise "offtarget"
    if fraction < 0.7:
        return 'translocation'
    else:
        return 'offtarget'
    

# Convert to DataFrame and save to CSV
translocations_df = pd.DataFrame(translocation_rows)
translocations_df.to_csv('translocations_debug.csv', index = False)
#differentiate translocations and off target amplifications
translocations_df['classification'] = translocations_df.apply(classify_translocation, axis=1)
#separate into two dataframes (translocations and offtargets)
df_translocation = translocations_df[translocations_df['classification'] == 'translocation']
df_offtarget = translocations_df[translocations_df['classification'] == 'offtarget']
df_translocation.to_csv('potential_translocations.csv', index=False)
df_offtarget.to_csv('potential_offtargets.csv', index=False)


# Create DataFrame
dfmain_mod = pd.DataFrame(rows)

# Rename the columns 'Unnamed: 1' and 'Unnamed: 2' to 'bc' and 'umi' in both DataFrames
dfmain_mod.rename(columns={'Unnamed: 0': 'bc', 'Unnamed: 1': 'umi'}, inplace=True)

# Save the DataFrame to CSV
dfmain_mod.to_csv('bc_umi_mod_seq.csv', index=False)

# Define a function to calculate HDR count
def calculate_hdr_count(hdr_bc_string):
    return sum(int(hdrbc.split('_')[1]) for hdrbc in hdr_bc_string.split(';'))

if os.path.exists(filename_hdr):
    # Apply the function to the 'hdr_bc' column to create a new 'HDR' column
    dfmain_hdr['HDR'] = dfmain_hdr['hdr_bc'].apply(calculate_hdr_count)


    # Write the updated DataFrame back to CSV without the index
    dfmain_hdr.to_csv('bc_umi_hdr_seq.csv', index=False)

    # Perform a full join (equivalent to outer join in pandas) on the 'bc' and 'umi' columns
    merged_df = pd.merge(dfmain_hdr, dfmain_mod, how='outer', on=['bc', 'umi'])
    if 'raw_bc_x' in merged_df.columns or 'raw_bc_y' in merged_df.columns:
        merged_df['raw_bc'] = merged_df.get('raw_bc_x', pd.Series(dtype=object)).combine_first(
            merged_df.get('raw_bc_y', pd.Series(dtype=object))
        )
        merged_df = merged_df.drop(columns=[column for column in ['raw_bc_x', 'raw_bc_y'] if column in merged_df.columns])

    # Write the merged DataFrame to a CSV file
    merged_df.to_csv("Hdr_mods.csv", index=False, na_rep="")
else:
    #add missing HDR columns to dataframe and save as Hdr_mods.csv
    dfmain_mod['HDR'] = ""
    dfmain_mod["hdr_bc"] = ""
    dfmain_mod.to_csv('Hdr_mods.csv', index=False)
