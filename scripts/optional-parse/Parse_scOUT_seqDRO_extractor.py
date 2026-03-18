# -*- coding: utf-8 -*-
"""
Assign detailed repair outcomes from parse data to mouse cells edited at B2m or Son
"""
#specify target as first command line argument

#change WT_amplicon and left_anchor/right_anchor/HDR_anchor according to targeted gene
#change sub_library according to barcode/sublibrary specification
#change hdr barcode length (hdrbc_len, 12/15) according to HDR template specification
#change zero_pad_barcode() depending if your barcodes are already zero padded
#don't forget to specify the read length

import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
from collections import defaultdict
import os
import csv
import re
import math
import sys
from pathlib import Path

DRO_ALLOWED_FILES = {
    "AlleleFrequencies.xlsx",
    "AlleleFrequencies_sized.xlsx",
    "EditingOutcomesByCellBarcode.xlsx",
    "RepairOutcomeTable.csv",
    "UMIRepairOutcomes.csv",
    "UMI_per_cell_histogram.png",
    "aln_score_filtered_out_data.csv",
    "EditingOutcomeAssignment.log",
    "gap_filtered_out_data.csv",
    "reads_per_cell_histogram.png",
    "AllelesByCellBarcode.csv",
    "AlleleCount.csv",
}


def cleanup_dro_outputs() -> None:
    for path in Path(".").iterdir():
        if path.is_file() and path.name not in DRO_ALLOWED_FILES:
            path.unlink()


target = sys.argv[1]
read_length = int(sys.argv[2])
cbc_path = sys.argv[3]
sub_library = ''


#define scOUT target (GAPDH_sg13_SNP, Son_sg4, B2m_sg1, GAPDH_sg13)

if target not in ["GAPDH_sg13_SNP", "Son_sg4", "B2m_sg1", "GAPDH_sg13"]:
    sys.exit("Unrecognized scOUT target. Aborting.")
elif target == "Son_sg4":
        # Define the anchors - Son_sg4
    left_anchor = 'CACCCAGC' #PAM proximal
    right_anchor = 'CGCTATTTTG' #PAM distal
    HDR_anchor = 'GGCATGC' #PAM proximal including PAM mutation if applicable
    #Son_sg4 150 bp
    WT_amplicon = 'GATAGACGTAAATAAAAATGCTGTAACCGACTTATCTAATAAAAATTGGCACCCAGCCGCTATTTTGTTGACTGAGGAAGTTTATGTTAATTTTTTAGGGTCTGATAGAATATTCATGTGTATTACAGTGGTATTCATATGCTATGTCTCT'
    hdrbc_len = 15
    """thresholds = {
    'D1_D40': 60,
    'D40_D50': 60,
    'D50_D60': 60,
    'D60plus': 40,
    'I1_I20': 60,
    'I20_I35': 45,
    'I35_I50': 35,
    'I50plus': 30
    }"""
elif target == "B2m_sg1":
    # Define the anchors - B2m_sg1 (other direction than son)
    left_anchor = 'ACTTGGAT' #PAM proximal
    right_anchor = 'ACTTCTCATT' #PAM distal
    HDR_anchor = 'TCAATGCAGT' #PAM proximal including PAM mutation if applicable
    #B2m_sg1 150 bp
    WT_amplicon = 'GTATTTTGATCAGAATAATAAATATAATTTTAAGAACAATAGTTGATCATATGCCAAACCCTCTGTACTTCTCATTACTTGGATGCAGTTACTCATCTTTGGTCTATCACAACATAAGTGACATACTTTCCTTTTGGTAAAGCAAAGAGGC'
    hdrbc_len = 15
    """thresholds = {
    'D1_D40': 60,
    'D40_D50': 60,
    'D50_D60': 60,
    'D60plus': 40,
    'I1_I20': 60,
    'I20_I35': 45,
    'I35_I50': 35,
    'I50plus': 30
    }"""
elif target == "GAPDH_sg13":
    # Define the anchors - GAPDH_sg13 (other direction than son)
    left_anchor = 'CTCCTCAC' #PAM proximal
    right_anchor = 'AGTTGCCATG' #PAM distal
    HDR_anchor = 'CTCCCCTTGT' #PAM proximal including PAM mutation
    #GAPDH sg13 151 bp
    WT_amplicon = 'GTCCCTGCCACACTCAGTCCCCCACCACACTGAATCTCCCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGTCATGTACCATCAATAAAGTACCCTGTGCTCAACCAGTTACTTGTCCTGT'
    hdrbc_len = 15
    """thresholds = {
    'D1_D40': 40,
    'D40_D50': 40,
    'D50_D60': 30,
    'D60plus': 24,
    'D90plus': 100, #ignore deletions over 90nt
    'I1_I20': 45,
    'I20_I35': 45,
    'I35_I50': 35,
    'I50plus': 30
    }"""
    if read_length == 99:
        WT_amplicon = 'GTCCCTGCCACACTCAGTCCCCCACCACACTGAATCTCCCCTCCTCACAGTTGCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCT'
        """thresholds = {
        'D1_D40': 60,
        'D40_D50': 40,
        'D50_D60': 40,
        'D60plus': 35,
        'D90plus': 100, #ignore deletions over 90nt
        'I1_I20': 44,
        'I20_I35': 40,
        'I35_I50': 35,
        'I50plus': 30
        }"""
elif target == "GAPDH_sg13_SNP":
    #GAPDH_sg13_SNP
    left_anchor = 'CTCCTCAC'  # PAM proximal
    right_anchor = 'AGTTTCCATG'  # PAM distal
    HDR_anchor = 'CTCCCCTACT'  # PAM proximal including PAM mutation
    #GAPDH SNP 102 bp
    WT_amplicon = 'GTCCCTGCCACACTCAGTCCCCCACCACACTGAATCTCCCCTCCTCACAGTTTCCATGTAGACCCCTTGAAGAGGGGAGGGGCCTAGGGAGCCGCACCTTGT'
    hdrbc_len = 12
    """thresholds = {
    'D1_D40': 60,
    'D40_D50': 40,
    'D50_D60': 40,
    'D60plus': 35,
    'D90plus': 100, #ignore deletions over 90nt
    'I1_I20': 44,
    'I20_I35': 40,
    'I35_I50': 35,
    'I50plus': 30
    }"""




def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def flatten(l):
    return [item for sublist in l for item in sublist]

ROTable = pd.read_csv("../Hdr_mods.csv", header=0, low_memory=False)

# Number of nuclei barcodes detected before filtering
print(f"Number of unique barcodes before filtering: {ROTable['bc'].nunique()}")
def reformat_barcodes(df, col_name, sub_library):
    def reformat(x):
        # OPTIONAL: Reformat the barcode by replacing colons with underscores
        #formatted_barcode = '_'.join([y.zfill(2) if y != '-' else '-' for y in x.split(':')])
        
        # Append sub_library if it's defined and non-empty
        formatted_barcode = x
        if sub_library:
            formatted_barcode += f"__{sub_library}"
        return formatted_barcode
    
    df[col_name] = df[col_name].apply(reformat)
    return df
if sub_library:
    ROTable = reformat_barcodes(ROTable, 'bc', sub_library)

    


# Read the cell metadata
cell_metadata = pd.read_csv(f"{cbc_path}cell_metadata_unfiltered.csv", header=0, low_memory=False)
whitelist_unfiltered = set(cell_metadata['bc_wells'].tolist())
cell_metadata.set_index('bc_wells', inplace=True)


all_barcodes_before_filtering = set(ROTable['bc'].unique())

whitelist_df = pd.read_csv(f"{cbc_path}cbcs_allsamples_filtered.csv", header=0, low_memory=False)
whitelist = set(whitelist_df['Cells.all_samples.data.'].tolist())

barcodes_filtered = set(whitelist)
filtered_overlap = all_barcodes_before_filtering.intersection(barcodes_filtered)
num_filtered_overlap = len(filtered_overlap)
print(f"Number of unique barcodes in filtered whitelist: {num_filtered_overlap}")

ROTable = ROTable.fillna(0)


# Define a class for Trie Node
class TrieNode:
    def __init__(self):
        self.children = {}
        self.is_end_of_word = False
        self.word = None

# Define a class for the Trie itself
class Trie:
    def __init__(self):
        self.root = TrieNode()

    # Insert a barcode into the Trie
    def insert(self, word):
        node = self.root
        parts = word.split('_')
        for part in parts:
            if part not in node.children:
                node.children[part] = TrieNode()
            node = node.children[part]
        node.is_end_of_word = True
        node.word = word

    # Search for matching barcodes in the Trie
    def search(self, barcode):
        node = self.root
        parts = barcode.split('_')
        return self._search_helper(node, parts, 0)

    # Helper function for search
    def _search_helper(self, node, parts, idx):
        if idx == len(parts):
            return [node.word] if node.is_end_of_word else []
        
        matches = []
        part = parts[idx]
        
        if part == '-':
            # Wildcard case: Check all children
            for child_part, child_node in node.children.items():
                matches += self._search_helper(child_node, parts, idx + 1)
        elif part in node.children:
            # Exact match case
            matches += self._search_helper(node.children[part], parts, idx + 1)
        
        return matches

# Preprocess whitelist into a Trie
def preprocess_whitelist_trie(whitelist):
    trie = Trie()
    for item in whitelist:
        trie.insert(item)
    return trie

# Optimized barcode matching function using Trie
def barcode_matches_optimized_trie(barcode, trie):
    matches = trie.search(barcode)
    if len(matches) == 1:
        print(f"Original barcode: {barcode}, Corrected barcode: {matches[0]}")
        return matches[0]
    return None

# Update DataFrame using Trie for barcode correction
def update_dataframe_barcodes_optimized_trie(df, col_name, whitelist):
    trie = preprocess_whitelist_trie(whitelist)
    
    # Extract all unique barcodes that contain '-'
    unique_barcodes_with_dash = df[col_name][df[col_name].str.contains('-')].unique()
    
    # Create a mapping dictionary
    barcode_mapping = {}
    for barcode in unique_barcodes_with_dash:
        corrected_barcode = barcode_matches_optimized_trie(barcode, trie)
        barcode_mapping[barcode] = corrected_barcode

    # Update the DataFrame using the mapping
    df[col_name] = df[col_name].map(lambda x: barcode_mapping.get(x, x))
    
    return df

# Now you can use the updated function with the Trie implementation
ROTable_updated = update_dataframe_barcodes_optimized_trie(ROTable, 'bc', whitelist_unfiltered)

#remove NA cells
ROTable_updated = ROTable_updated[ROTable_updated['bc'].notna()]




# Get the count of unique values in these filtered rows
rows_with_hyphen = ROTable_updated[ROTable_updated['bc'].str.contains('-')]
unique_count = rows_with_hyphen['bc'].nunique()
print("Number of unique rows containing incomplete barcode:", unique_count)

ROTable = ROTable_updated

# Consolidate identical UMIs (need to adjust this part)
def merge_strings_mods(s):
    counts = defaultdict(lambda: [0, 0])  # This will hold two counts per key
    for item in s.dropna():
        parts = str(item).split(';')
        for part in parts:
            if part != '0':
                try:
                    mod, count, seq, seq_count, alnrefseq = part.split('_')
                    count = int(count)
                    seq_count = int(seq_count)
                    counts[(mod, seq)][0] += count
                    counts[(mod, seq)][1] += seq_count
                except ValueError:
                    print(f"Skipping part with unexpected format: {part}")

    return ';'.join(f'{mod}_{count_vals[0]}_{seq}_{count_vals[1]}_{alnrefseq}' for (mod, seq), count_vals in counts.items())

def merge_strings_hdrbc(s):
    counts = defaultdict(lambda: [0, 0])  # This will hold two counts per key
    for item in s.dropna():
        parts = str(item).split(';')
        for part in parts:
            if part != '0':
                try:
                    mod, count, seq, seq_count, hd = part.split('_')
                    count = int(count)
                    seq_count = int(seq_count)
                    counts[(mod, seq)][0] += count
                    counts[(mod, seq)][1] += seq_count
                except ValueError:
                    print(f"Skipping part with unexpected format: {part}")

    return ';'.join(f'{mod}_{count_vals[0]}_{seq}_{count_vals[1]}_{hd}' for (mod, seq), count_vals in counts.items())

def merge_strings_wtseq(s):
    counts = defaultdict(lambda: [0, 0])  # This will hold two counts per key
    for item in s.dropna():
        parts = str(item).split(';')
        for part in parts:
            if part != '0':
                try:
                    seq, seq_count = part.split('_')
                    seq_count = int(seq_count)
                    counts[seq][0] += seq_count
                except ValueError:
                    print(f"Skipping part with unexpected format: {part}")

    return ';'.join(f'{seq}_{count_vals[0]}' for seq, count_vals in counts.items())


def sum_numeric(s):
    return s.sum()
ROTable = ROTable.groupby(['bc', 'umi']).agg({
    'hdr_bc': merge_strings_hdrbc,
    'HDR': 'sum',
    'Reference_UNMODIFIED': 'sum',
    'Reference_MODIFIED': 'sum',
    'mods': merge_strings_mods,
    'wt_seq': merge_strings_wtseq,
    'avg_aln_score': 'mean'
    # Add other columns as needed
}).reset_index()


#assign UMI_repair Outcome
conditions = [
    (ROTable['Reference_UNMODIFIED'] >= 1) & (ROTable['Reference_UNMODIFIED'] > ROTable['HDR']) & (ROTable['Reference_UNMODIFIED'] > ROTable['Reference_MODIFIED']),
    (ROTable['HDR'] >= 1) & (ROTable['HDR'] > ROTable['Reference_UNMODIFIED']) & (ROTable['HDR'] > ROTable['Reference_MODIFIED']),
    (ROTable['Reference_MODIFIED'] >= 1) & (ROTable['Reference_MODIFIED'] > ROTable['Reference_UNMODIFIED']) & (ROTable['Reference_MODIFIED'] > ROTable['HDR']),
    ]
values = ['WT', 'HDR', 'NHEJ']
ROTable['UMI_RepairOutcome'] = np.select(conditions, values, default='NA')
#filter NA UMIs
ROTable = ROTable[ROTable['UMI_RepairOutcome'] != 'NA']
ROTable = ROTable.fillna(0)

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

def should_filter(row):
    mods = str(row['mods']).split(';')[0]
    # Take only the first modification entry
    first_mod = mods.split(';')[0]
    # Find all continuous gaps (sequences of dashes)
    gaps = first_mod.split('_')
    if len(gaps) >= 2:
        mod = gaps[0]
        seq = gaps[2]
        #continuous_gap_lengths = count_continuous_gaps(seq)

        if row['UMI_RepairOutcome'] == 'NHEJ':
            return all_gaps_trailing(gaps[2], mod)
            #return (check_gaps(continuous_gap_lengths) or all_gaps_trailing(gaps[2], mod))
        else:
            return False
    else:
        return False


mask = ROTable.apply(should_filter, axis=1)
filtered_out_df = ROTable[mask]
filtered_out_df.to_csv('gap_filtered_out_data.csv', index=False)

# Apply the filter function to each row
ROTable = ROTable[~ROTable.apply(should_filter, axis=1)]


#filter NHEJ UMIs without mods
condition = (ROTable['UMI_RepairOutcome'] == 'NHEJ') & (ROTable['mods'] == 0)
ROTable = ROTable[~condition]


#filter to only cell barcodes that occur in the RNA-seq - OPTIONAL
#ROTable = ROTable[ROTable['bc'].isin(whitelist)]


#filter NHEJ UMIs without clear mod (most likely bad alignment) and with bad alignment score
def check_reads(row):
    # Check if UMI_RepairOutcome is 'NHEJ'
    if row['UMI_RepairOutcome'] == 'WT':
        parts = str(row['wt_seq']).split(';')[0]
        aln_score = float(row['avg_aln_score'])
        min_required = 0.4 * row['Reference_UNMODIFIED']
        try:
            number = int(parts.split('_')[1])
            if number < min_required:
                return False
            if aln_score < 55:
                return False
            return True
        except (IndexError, ValueError):
            return True
    elif row['UMI_RepairOutcome'] == 'NHEJ':
        # Splitting the sequences by semicolon
        parts = str(row['mods']).split(';')[0]
        aln_score = float(row['avg_aln_score'])
        # Getting the minimum required reads (70% of the reads column)
        min_required = 0.6 * row['Reference_MODIFIED']
        # Extract the number after the first underscore and before the second underscore
        # Safe guard to handle any unexpected parts without the expected underscore format
        try:
            number = int(parts.split('_')[1])
            mod = parts.split('_')[0]
            
            if 'S' in mod and (mod.startswith('D') or mod.startswith('I')):
                mod = mod.split('S')[0]  # Strip out the 'S' and the following substitution size

            
            if number < min_required:
                return False
            return True
        except (IndexError, ValueError):
            return True
        return True
    else:
        return True



mask = ROTable.apply(check_reads, axis=1)
filtered_out_df = ROTable[~mask]
filtered_out_df.to_csv('aln_score_filtered_out_data.csv', index=False)

# Applying the function to filter the DataFrame
ROTable = ROTable[ROTable.apply(check_reads, axis=1)]


#filter to only cell barcodes that occur in the RNA-seq - OPTIONAL
#ROTable = ROTable[ROTable['bc'].isin(whitelist)]


#filter out NHEJ UMIs without mods
ROTable = ROTable[~((ROTable['UMI_RepairOutcome'] == 'NHEJ') & (ROTable['mods'].isna() | (ROTable['mods'] == '')))]

ROTable.to_csv('RepairOutcomeTable.csv', index=False)




ROTable.groupby('bc')['UMI_RepairOutcome'].value_counts().to_csv('UMIRepairOutcomes.csv')
UMIRepairOutcome = ROTable.groupby('bc')['UMI_RepairOutcome'].value_counts().to_dict()

#UMI per Cell histogram
UMIperBC = ROTable['bc'].value_counts()
plt.hist(UMIperBC, bins=range(0, int(max(UMIperBC)) + 2, 2), edgecolor='black')
plt.xlim(0, max(UMIperBC)+2)
# Setting labels and title
plt.xlabel('UMI per Barcode')
plt.ylabel('Frequency')
plt.title('Histogram of UMI per Barcode')
# Saving the plot
plt.savefig('UMI_per_cell_histogram.png', format='png', dpi=300)
plt.close()


print('median UMIs per cell:', ROTable['bc'].value_counts().median())
# Calculate the sum across the columns HDR, Reference_UNMODIFIED, and Reference_MODIFIED
ROTable['totalrc'] = ROTable[['HDR', 'Reference_UNMODIFIED', 'Reference_MODIFIED']].sum(axis=1)


#Reads per Cell histogram
reads_per_bc = ROTable.groupby('bc')['totalrc'].sum()
plt.hist(reads_per_bc, bins=range(0, int(max(reads_per_bc)) + 100, 100), edgecolor='black')
plt.xlabel('Reads per Barcode')
plt.ylabel('Frequency')
plt.title('Histogram of Reads per Barcode')
plt.savefig('reads_per_cell_histogram.png', format='png', dpi=300)
plt.close()

#Reads per UMI histogram
reads_per_UMI = ROTable['totalrc']
plt.hist(reads_per_UMI, bins=range(0, int(max(reads_per_UMI)) + 100, 100), edgecolor='black')
plt.xlabel('Reads per UMI')
plt.ylabel('Frequency')
plt.title('Histogram of Reads per UMI')
plt.savefig('reads_per_UMI_histogram.png', format='png', dpi=300)
plt.close()

# Group the table by the 'bc' column and calculate the mean of the 'sum_cols' for each 'bc'
meanrc_per_cell = ROTable.groupby('bc')['totalrc'].mean()
print('mean read count per cell:', meanrc_per_cell.mean())




#filter out true modifications
UMImod = dict()
UMIhdrbc = dict()
UMIwt = dict()

for index, row in ROTable.iterrows():
    if type(row['mods']) == str and row['mods'] != '0' and row['UMI_RepairOutcome'] == 'NHEJ':
        Mods=str(row['mods']).split(sep=';')
        if int(row['Reference_MODIFIED']) > int(row['HDR']):
             ModsReadCount=int(row['Reference_MODIFIED'])
        elif int(row['Reference_MODIFIED']) < int(row['HDR']):
             ModsReadCount=int(row['HDR'])
        TotalCount_Mod=0
        TotalCount_Modseq=0
        ModsDict={}
        TrueMods=list()
        for Mod in Mods:
            Mod, count, modseq, modseq_count, refalnseq = Mod.split(sep='_')
            ModsDict[Mod] = {'count' : int(count), 'modseq' : modseq, 'modseq_count' : modseq_count, 'refalnseq' : refalnseq}
            TotalCount_Mod += int(count)
            TotalCount_Modseq += int(modseq_count)
            if int(count)/ModsReadCount > 0.5 and int(count) >= 1:
                TrueMods.append([Mod, int(count), ModsReadCount, modseq, int(modseq_count), refalnseq])
        UMImod[row['bc'], row['umi']] = TrueMods
    elif type(row['hdr_bc']) == str and row['hdr_bc'] != '0' and row['UMI_RepairOutcome'] == 'HDR':
        HdrBarcodes=str(row['hdr_bc']).split(sep=';')
        HdrReadCount=int(row['HDR'])
        HdrBarcodeDict={}
        TrueHdrBarcodes=list()
        for HdrBc in HdrBarcodes:
            HdrBarcode, count, hdrbc_seq, hdrbc_seq_count, hdrbc_hd = HdrBc.split(sep='_')
            HdrBarcodeDict[HdrBarcode] = {'count' : int(count), 'hdrbc_seq' : hdrbc_seq, 'hdrbc_seq_count' : hdrbc_seq_count, 'hdrbc_hd' : hdrbc_hd}
            if int(count)/HdrReadCount > 0.5 and int(count) >= 1:
                TrueHdrBarcodes.append([HdrBarcode, int(count), HdrReadCount, hdrbc_seq, int(hdrbc_seq_count), hdrbc_hd])
        UMIhdrbc[row['bc'], row['umi']] = TrueHdrBarcodes
    elif row['UMI_RepairOutcome'] == 'WT':
        wt_seq, wt_seq_count =row['wt_seq'].split("_") #use wt seq count?
        UMIwt[row['bc'], row['umi']] = ['WT', int(row['Reference_UNMODIFIED']), wt_seq]
    else:
        print(f"No valid UMI repair outcome: {row['bc']}, {row['umi']}")
        continue


#count UMIs AND reads per hdrbc and add value
Cell_hdrbc = dict()
for key, value in UMIhdrbc.items():
    if key[0] not in Cell_hdrbc.keys() and len(value) != 0:
        hdrbc = value[0][0]
        hdrbc_seq = value[0][3]
        hdrbc_umicount = 1
        hdrrc = value[0][2]
        hdrbc_hd = value[0][5]
        hdrbc_dict = {'umi_count' : hdrbc_umicount, 'rc' : hdrrc, 'seq' : hdrbc_seq, 'hd' : hdrbc_hd}
        Cell_hdrbc[key[0]] = {'hdrbcs' : {hdrbc : hdrbc_dict}}
    elif key[0] in Cell_hdrbc.keys() and len(value) != 0:
        hdrbc = value[0][0]
        hdrbc_seq = value[0][3]
        hdrrc = value[0][2]
        if hdrbc in Cell_hdrbc[key[0]]['hdrbcs'].keys():
            hdrbc_umicount = Cell_hdrbc[key[0]]['hdrbcs'][hdrbc]['umi_count'] + 1
            hdrrc = Cell_hdrbc[key[0]]['hdrbcs'][hdrbc]['rc'] + hdrrc
            Cell_hdrbc[key[0]]['hdrbcs'][hdrbc]['umi_count'] = hdrbc_umicount
            Cell_hdrbc[key[0]]['hdrbcs'][hdrbc]['rc'] = hdrrc
        elif hdrbc not in Cell_hdrbc[key[0]]['hdrbcs'].keys():
            hdrbc_umicount = 1
            hdrbc_dict = {'umi_count' : hdrbc_umicount, 'rc' : hdrrc, 'seq' : hdrbc_seq, 'hd' : hdrbc_hd}
            Cell_hdrbc[key[0]]['hdrbcs'][hdrbc] = hdrbc_dict
        else:
            print("This UMI is in line 330:", value)
    else:
        print("This UMI is in line 332:", value)
 
    

#create dict with total readcounts per cell
Cell_totalrc = ROTable.groupby('bc')['totalrc'].sum().to_dict()
#create dict with total UMIs per cell
Cell_totalumi = ROTable.groupby('bc').size().to_dict()

#first remove all barcodes, then add barcode with same information and new umi count
#add HDR readcount here if necessary!
hdrbcs_to_remove = {}
for key, value in Cell_hdrbc.items():
    if len(value['hdrbcs']) == 2:
        if hamming_distance(list(value['hdrbcs'].keys())[0], list(value['hdrbcs'].keys())[1]) == 1:
            min_hdrbc = min(value['hdrbcs'].items(), key=lambda x: x[1]['umi_count'])
            max_hdrbc = max(value['hdrbcs'].items(), key=lambda x: x[1]['umi_count'])
            new_umicount = sum(item[1]['umi_count'] for item in value['hdrbcs'].items())
            max_hdrbc[1]['umi_count'] = new_umicount
            hdrbcs_to_remove[key] = max_hdrbc
for bc, hdrbc in hdrbcs_to_remove.items():
    del Cell_hdrbc[bc]['hdrbcs']
    Cell_hdrbc[bc]['hdrbcs']= {hdrbc[0]:hdrbc[1]}
 
#count UMIs AND reads per mod and add value
Cell_mod = dict()
for key, value in UMImod.items():
    if len(value) > 1:
        modseqs = set()
        modalnrefseqs = set()
        mods = list()
        for l in value:
            modseqs.add(l[3])
            modalnrefseqs.add(l[5])
            mods.append(l[0])
        mods = ''.join(sorted(mods))
        if len(modseqs) == 1:
            if key[0] not in Cell_mod.keys():
                mod = mods
                mod_seq = str(list(modseqs)[0])
                mod_alnrefseq = str(list(modalnrefseqs)[0])
                mod_count = 1
                mod_rc = value[0][2]
                mod_dict = {'umi_count' : mod_count, 'rc' : mod_rc, 'seq' : mod_seq, 'alnrefseq' : mod_alnrefseq}
                Cell_mod[key[0]] = {'mods' : {mod : mod_dict}}
            elif key[0] in Cell_mod.keys():
                mod = mods
                mod_seq = str(list(modseqs)[0])
                mod_alnrefseq = str(list(modalnrefseqs)[0])
                mod_rc = value[0][2]
                if mod in Cell_mod[key[0]]['mods'].keys():
                    mod_count = Cell_mod[key[0]]['mods'][mod]['umi_count'] + 1
                    Cell_mod[key[0]]['mods'][mod]['umi_count'] = mod_count
                    mod_rc = Cell_mod[key[0]]['mods'][mod]['rc'] + mod_rc
                    Cell_mod[key[0]]['mods'][mod]['rc'] = mod_rc
                elif mod not in Cell_mod[key[0]]['mods'].keys():
                    mod_count = 1
                    mod_dict = {'umi_count' : mod_count, 'rc' : mod_rc, 'seq' : mod_seq, 'alnrefseq' : mod_alnrefseq}
                    Cell_mod[key[0]]['mods'][mods] = mod_dict
                else:
                    print("This UMI is in line 388", value)
        else:
            for l in value:
                if key[0] not in Cell_mod.keys():
                    mod = l[0]
                    mod_seq = str(l[3])
                    mod_alnrefseq = str(l[5])
                    mod_count = 1
                    mod_rc = value[0][2]
                    mod_dict = {'umi_count' : mod_count, 'rc' : mod_rc, 'seq' : mod_seq, 'alnrefseq' : mod_alnrefseq}
                    Cell_mod[key[0]] = {'mods' : {mod : mod_dict}}
                elif key[0] in Cell_mod.keys():
                    mod = l[0]
                    mod_seq = str(l[3])
                    mod_alnrefseq = str(l[4])
                    mod_rc = l[2]
                    if mod in Cell_mod[key[0]]['mods'].keys():
                        mod_count = Cell_mod[key[0]]['mods'][mod]['umi_count'] + 1
                        mod_rc = Cell_mod[key[0]]['mods'][mod]['rc'] + mod_rc
                        Cell_mod[key[0]]['mods'][mod]['umi_count'] = mod_count
                        Cell_mod[key[0]]['mods'][mod]['rc'] = mod_rc
                    elif mod not in Cell_mod[key[0]]['mods'].keys():
                        mod_count = 1
                        mod_rc = l[2]
                        mod_dict = {'umi_count' : mod_count, 'rc' : mod_rc, 'seq' : mod_seq, 'alnrefseq' : mod_alnrefseq}
                        Cell_mod[key[0]]['mods'][mod] = mod_dict
                    else:
                        print("This UMI is in line 413", value)
                else:
                    print("This UMI is in line 415", value)
    elif len(value) == 1:
        if key[0] not in Cell_mod.keys():
            mod = value[0][0]
            mod_seq = str(value[0][3])
            mod_alnrefseq = str(value[0][5])
            mod_count = 1
            mod_rc = value[0][2]
            mod_dict = {'umi_count' : mod_count, 'rc' : mod_rc, 'seq' : mod_seq, 'alnrefseq' : mod_alnrefseq}
            Cell_mod[key[0]] = {'mods' : {mod : mod_dict}}
        elif key[0] in Cell_mod.keys():
            mod = value[0][0]
            mod_seq = str(value[0][3])
            mod_alnrefseq = str(value[0][5])
            mod_rc = value[0][2]
            if mod in Cell_mod[key[0]]['mods'].keys():
                mod_count = Cell_mod[key[0]]['mods'][mod]['umi_count'] + 1
                mod_rc = Cell_mod[key[0]]['mods'][mod]['rc'] + mod_rc
                Cell_mod[key[0]]['mods'][mod]['umi_count'] = mod_count
                Cell_mod[key[0]]['mods'][mod]['rc'] = mod_rc
            elif mod not in Cell_mod[key[0]]['mods'].keys():
                mod_count = 1
                mod_rc = value[0][2]
                mod_dict = {'umi_count' : mod_count, 'rc' : mod_rc, 'seq' : mod_seq, 'alnrefseq' : mod_alnrefseq}
                Cell_mod[key[0]]['mods'][mod] = mod_dict
            else:
                print("This UMI is in line 439", value)
        else:
            print("This UMI is in line 441", value)
    elif len(value) == 0:
        continue
    else:
        print("This UMI is in line 445", value)

#count UMIs per wt and rc per wt and add value
Cell_wt = dict()
for key, value in UMIwt.items():
    if value[1] >= 1:
        if key[0] not in Cell_wt.keys():
            wt_seq = value[2]
            wt_count = 1
            wt_rc = value[1]
            wt_dict = {'umi_count' : wt_count, 'rc' : wt_rc, 'seq' : wt_seq}
            Cell_wt[key[0]] = wt_dict
        elif key[0] in Cell_wt.keys():
            wt_seq = value[2]
            wt_count = Cell_wt[key[0]]['umi_count'] + 1
            wt_rc = value[1]
            wt_rc = Cell_wt[key[0]]['rc'] + wt_rc
            Cell_wt[key[0]]['umi_count'] = wt_count
            Cell_wt[key[0]]['rc'] = wt_rc
        else:
            print("This UMI is in line ", value)



#output Cell_hdrbc_mod only with important things: rc_perc, rc, modseq ==> for each allele
Cell_hdrbc_mod = {}
NHEJperclist = []


hdrbc_perc_list = []
mod_perc_list = []
wt_perc_list = []

totalUMIcounts = dict()


for bc in ROTable['bc'].unique():
    nhej_rc_perc = 0
    hdr_rc_perc = 0
    wt_rc_perc = 0
    cell_rc = Cell_totalrc[bc]
    Cell_hdrbc_mod[bc] = {'alleles' : dict()}
    if bc in Cell_hdrbc.keys():
        UMI_counts_hdrbc = {key: value['umi_count'] for key, value in Cell_hdrbc[bc]['hdrbcs'].items()}
        hdr_rc_dict = {key: value['rc'] for key, value in Cell_hdrbc[bc]['hdrbcs'].items()}
        for hdrbc, hdr_rc in hdr_rc_dict.items():
            hdr_rc_perc = hdr_rc / cell_rc * 100
            hdrbc_perc_list.append(hdr_rc_perc)
            Cell_hdrbc[bc]['hdrbcs'][hdrbc].update({'rc_perc' : hdr_rc_perc})
            Cell_hdrbc_mod[bc]['alleles'].update(Cell_hdrbc[bc]['hdrbcs'])
            #OPTIONAL:filter by read count percentage
            """if hdr_rc_perc > 20:
                Cell_hdrbc[bc]['hdrbcs'][hdrbc].update({'rc_perc' : hdr_rc_perc})
                Cell_hdrbc_mod[bc]['alleles'].update(Cell_hdrbc[bc]['hdrbcs'])
            else:
                print("hdr_rc_perc < 20", hdr_rc_perc)"""
    if bc in Cell_mod.keys():
        UMI_counts_mods = {key: value['umi_count'] for key, value in Cell_mod[bc]['mods'].items()}
        mod_rc_dict = {key: value['rc'] for key, value in Cell_mod[bc]['mods'].items()}
        for mod, mod_rc in mod_rc_dict.items():
            mod_rc_perc = mod_rc / cell_rc * 100
            mod_perc_list.append(mod_rc_perc)
            Cell_mod[bc]['mods'][mod].update({'rc_perc' : mod_rc_perc})
            Cell_hdrbc_mod[bc]['alleles'].update(Cell_mod[bc]['mods'])
            #OPTIONAL:filter by read count percentage
            """if mod_rc_perc > 10:
                Cell_mod[bc]['mods'][mod].update({'rc_perc' : mod_rc_perc})
                Cell_hdrbc_mod[bc]['alleles'].update(Cell_mod[bc]['mods'])
            else:
                print("mod_rc_perc < 10", mod_rc_perc)"""
    if bc in Cell_wt.keys():
        UMI_counts_wt = Cell_wt[bc]['umi_count']
        wt_rc = Cell_wt[bc]['rc']
        wt_rc_perc = wt_rc / cell_rc * 100
        wt_perc_list.append(wt_rc_perc)
        Cell_wt[bc].update({'rc_perc' : wt_rc_perc})
        Cell_hdrbc_mod[bc]['alleles'].update({'WT' : Cell_wt[bc]})
        #OPTIONAL:filter by read count percentage
        """if wt_rc_perc > 20:
            Cell_wt[bc].update({'rc_perc' : wt_rc_perc})
            Cell_hdrbc_mod[bc]['alleles'].update({'WT' : Cell_wt[bc]})
        else:
            print("wt_rc_perc < 20", wt_rc_perc)"""
    else:
       print("This bc is in ROTable but not in a repair outcome dict", bc)
       pass


allele_dist = []
for bc, value in Cell_hdrbc_mod.items():
    for alleles in value.values():
        allele_dist.append(len(alleles))

counter = Counter(allele_dist)
#plot allele number histogram
plt.close()
plt.xlim(0, max(allele_dist)+1)
plt.hist(allele_dist, bins = range(0, int(max(allele_dist)) + 1, 1))
plt.xlabel('Number of Alleles per Cell')
plt.ylabel('Frequency')
plt.title('Histogram of Allele Numbers per Cell')
plt.savefig('allele_count_histogram.png', format='png', dpi=300)
plt.close()

# calculate number of alleles assigned (not including intermediates)
allele_counts = {}
for key, value in Cell_hdrbc_mod.items():
    allele_value = value.get('alleles')
    if isinstance(allele_value, list):
        allele_counts[key] = len(allele_value)
    elif isinstance(allele_value, dict):
        allele_counts[key] = len(allele_value)
    else:
       print("This cell is in line 620", key)

#number of alleles histogram
plt.hist(allele_counts.values(), bins = range(0, int(max(allele_counts.values())) + 1, 1))
plt.xlabel('Number of Alleles per cell')
plt.ylabel('Frequency')
plt.title('Histogram of number of Alleles')
plt.savefig('allele_count_histogram_2.png', format='png', dpi=300)
plt.close()

print("Allele numbers: ", Counter(allele_counts.values()))

#read percentage histogram
plt.hist([hdrbc_perc_list, mod_perc_list, wt_perc_list], bins=10, label=['HDR', 'NHEJ-Mod', 'WT'])
plt.xlabel('Percentage of reads per Cell')
plt.ylabel('Frequency')
plt.title('Histogram of Read Percentage per Cell')
plt.legend()
plt.xlim(0, 100) # Set the x-axis limits
plt.savefig('read_perc_histogram.png', format='png', dpi=300)
plt.close()


#read percentage histogram - only NHEJ
plt.hist(mod_perc_list, bins=50)
plt.xlabel('Percentage of NHEJ reads per Cell')
plt.ylabel('Frequency')
plt.title('Histogram of NHEJ')
plt.savefig('read_perc_NHEJ_histogram.png', format='png', dpi=300)
plt.close()

#read percentage histogram - only wt
plt.hist(wt_perc_list, bins=50)
plt.xlabel('Percentage of WT reads per Cell')
plt.ylabel('Frequency')
plt.title('Histogram of WT')
plt.savefig('read_perc_WT_histogram.png', format='png', dpi=300)
plt.close()

#read percentage histogram - only hdr
plt.hist(hdrbc_perc_list, bins=50)
plt.xlabel('Percentage of HDR reads per Cell')
plt.ylabel('Frequency')
plt.title('Histogram of HDR')
plt.savefig('read_perc_HDR_histogram.png', format='png', dpi=300)
plt.close()

#remove empty dictionaries from Cell_hdrbc_mod[bc]['alleles']
Cell_hdrbc_mod = {bc: inner_dict for bc, inner_dict in Cell_hdrbc_mod.items() if inner_dict['alleles']}


missing_rc_perc_entries = []
cell_rc = Cell_totalrc[bc]
for bc, bc_data in Cell_hdrbc_mod.items():
    for allele, allele_data in bc_data['alleles'].items():
        if 'rc_perc' not in allele_data:
            # Add 'rc_perc' with the value of 0.1
            missing_rc_perc_entries.append((bc, allele))
            mod_rc = Cell_hdrbc_mod[bc]['alleles'][allele]['rc']
            cell_rc = Cell_totalrc[bc]
            mod_rc_perc = mod_rc / cell_rc * 100
            Cell_hdrbc_mod[bc]['alleles'][allele]['rc_perc'] = mod_rc_perc



DetailedRepairOutcome = {}

with open('ModSelection.csv', 'w', newline='') as f:
    # Initialize the CSV writer
    writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
    
    # Write the header dynamically
    dynamic_header = ["bc", "allele_count", "total_rc", "total_UMIcount"]
    max_alleles = max(len(Cell_hdrbc_mod[bc]['alleles']) for bc in Cell_hdrbc_mod)  # Get the maximum number of alleles
    
    # Add dynamic columns to the header for alleles
    for i in range(1, max_alleles + 1):
        dynamic_header.extend([
            f"Allele{i}", f"rc_{i}", f"rc_perc_{i}", f"UMIcount_{i}", f"Sequence_{i}", 
            f"HDist_HDRBC_{i}", f"AlnRefSeq_{i}"  # Adding AlnRefSeq for each allele
        ])
    # Add static columns to the header
    dynamic_header.extend(["Ploidy", "DetailedRepairOutcome", "Sample", "In_Whitelist"])
    writer.writerow(dynamic_header)
    
    # Iterate over the data
    for bc, alleles in Cell_hdrbc_mod.items():
        allele_keys = sorted(Cell_hdrbc_mod[bc]['alleles'])
        allele_count = len(allele_keys)
        sample = cell_metadata.loc[bc, 'sample'] if bc in cell_metadata.index else 'Unknown'
        in_whitelist = bc in whitelist
        
        # Gather the data for this row
        row = [bc, allele_count]
        
        # Calculate total_rc by summing all rc values
        total_rc = sum(Cell_hdrbc_mod[bc]['alleles'][allele]['rc'] for allele in allele_keys)
        row.append(total_rc)
        total_umicount = sum(Cell_hdrbc_mod[bc]['alleles'][allele]['umi_count'] for allele in allele_keys)
        row.append(total_umicount)
        # Add dynamic data for each allele
        for allele in allele_keys:
            row.extend([
                allele,
                Cell_hdrbc_mod[bc]['alleles'][allele]['rc'],
                Cell_hdrbc_mod[bc]['alleles'][allele]['rc_perc'],
                Cell_hdrbc_mod[bc]['alleles'][allele]['umi_count'],
                Cell_hdrbc_mod[bc]['alleles'][allele]['seq'],
                Cell_hdrbc_mod[bc]['alleles'][allele].get('hd', ""),
                Cell_hdrbc_mod[bc]['alleles'][allele].get('alnrefseq', "") 
            ])
        
        # Fill in missing columns for rows with fewer alleles than the maximum
        for _ in range(max_alleles - len(allele_keys)):
            row.extend([None, None, None, None, None, None, None])
        
        # Add static data to the row
        if allele_count == 1:
            if total_umicount == 1:
                ploidy = "haploid"
                detailed_repair_outcome = "{}_NA".format(allele_keys[0])
            else:
                ploidy = "diploid"
                detailed_repair_outcome = "{}_{}".format(allele_keys[0], allele_keys[0])
        elif allele_count == 2:
            ploidy = "diploid"
            detailed_repair_outcome = "{}_{}".format(allele_keys[0], allele_keys[1])
        elif allele_count == 3:
            ploidy = "triploid"
            detailed_repair_outcome = "{}_{}_{}".format(allele_keys[0], allele_keys[1], allele_keys[2])
        else:
            ploidy = "polyploid"
            detailed_repair_outcome = "_".join(allele_keys)
        
        # Add the detailed repair outcome to the dictionary
        DetailedRepairOutcome[bc] = detailed_repair_outcome.split('_')
        
        # Add remaining static columns
        row.extend([ploidy, detailed_repair_outcome, sample, in_whitelist])
        
        # Write the row to the CSV
        writer.writerow(row)


      
# Step 1: Create a dictionary with HDR, important to set hdr barcode length!
DetailedRepairOutcome_HDR = {}
for key, value in DetailedRepairOutcome.items():
    # Convert all HDR-bc character alleles to "HDR" or f'IncompleteHDR-{hd}'
    transformed_values = []
    for item in value:
        if isinstance(item, str) and len(item) == hdrbc_len:
            hd_value = int(Cell_hdrbc_mod[key]['alleles'][item].get('hd', 0))  # Safely get 'hd' and convert it to int
            if hd_value >= hdrbc_len - 1:
                transformed_values.append("HDR")
            else:
                transformed_values.append(f"IncompleteHDR-{hd_value}")
        else:
            transformed_values.append(item)
    
    DetailedRepairOutcome_HDR[key] = transformed_values

# Step 2: Flatten the list of all alleles across all cells
AllAlleles = []
for alleles in DetailedRepairOutcome_HDR.values():
    for allele in alleles:
        AllAlleles.append(allele)
            
# Step 3: Count occurrences of each allele
AlleleCount = Counter(AllAlleles)

# Step 4: Write the allele counts to a CSV file
with open('AlleleCount.csv', mode='w') as fp:
    fp.write('Allele,Frequency\n')
    for tag, count in AlleleCount.items():
        fp.write('{},{}\n'.format(tag, count))

# Optional: Print the allele counts for debugging
print(AlleleCount)


AlleleDict = {}
for outcome in set(AllAlleles):
    AlleleDict[outcome] = set()
    
    
for cell, outcomes in DetailedRepairOutcome_HDR.items():
    for outcome in outcomes:
        if outcome in AlleleDict.keys():
            AlleleDict[outcome].add(cell)
        else:
           print("This cell is in line 787", cell, outcomes)

pd.DataFrame.from_dict(AlleleDict, orient='index').to_csv('AllelesByCellBarcode.csv')


#count occurences of WT or HDR alleles with other alleles
HDR_X = []
WT_X = []
for cell, outcomes in DetailedRepairOutcome_HDR.items():
    if 'HDR' in outcomes:
        for outcome in outcomes:
            if outcome != 'HDR':
                HDR_X.append(outcome)
for cell, outcomes in DetailedRepairOutcome_HDR.items():
    if 'WT' in outcomes:
        for outcome in outcomes:
            if outcome != 'WT':
                WT_X.append(outcome)
                        
HDR_X_allelecount = Counter(HDR_X)
with open('HDR_X_allelecount.csv', mode='w') as fp:
    fp.write('Allele,freq\n')
    for tag, count in HDR_X_allelecount.items():
        fp.write('{},{}\n'.format(tag, count))
        
WT_X_allelecount = Counter(WT_X)
with open('WT_X_allelecount.csv', mode='w') as fp:
    fp.write('Allele,freq\n')
    for tag, count in WT_X_allelecount.items():
        fp.write('{},{}\n'.format(tag, count))


#count occurences of any NHEJ alleles with other alleles
NHEJ_X = []

for cell, outcomes in DetailedRepairOutcome_HDR.items():
    if len(outcomes) == 2:
        if outcomes[0] != 'HDR' and outcomes[0] != 'WT' and outcomes[0] != 'ND':
            NHEJ_X.append(outcomes[1])
        elif outcomes[1] != 'HDR' and outcomes[1] != 'WT' and outcomes[1] != 'ND':
            NHEJ_X.append(outcomes[0])
                        
NHEJ_X_allelecount = Counter(NHEJ_X)
with open('NHEJ_X_allelecount.csv', mode='w') as fp:
    fp.write('Allele,freq\n')
    for tag, count in NHEJ_X_allelecount.items():
        fp.write('{},{}\n'.format(tag, count))


# Convert to DataFrame
df = pd.read_csv('ModSelection.csv', low_memory=False)



def sort_alleles_by_umi(df):
    max_alleles = df['allele_count'].max()

    def sort_row(row):
        allele_data = [
            {
                'Allele': row[f'Allele{i}'],
                'rc': row[f'rc_{i}'],
                'rc_perc': row[f'rc_perc_{i}'],
                'Sequence': row[f'Sequence_{i}'],
                'UMIcount': row[f'UMIcount_{i}'],
                'HDist_HDRBC': row[f'HDist_HDRBC_{i}'],
                'AlnRefSeq': row[f'AlnRefSeq_{i}']
            }
            for i in range(1, max_alleles + 1)
            if not pd.isna(row[f'UMIcount_{i}'])
        ]

        # Sort the groups by UMIcount in descending order
        allele_data.sort(key=lambda x: x['UMIcount'], reverse=True)

        # Create a sorted row
        sorted_row = {}
        for i in range(1, max_alleles + 1):
            if i <= len(allele_data):
                sorted_row[f'Allele{i}'] = allele_data[i-1]['Allele']
                sorted_row[f'rc_{i}'] = allele_data[i-1]['rc']
                sorted_row[f'rc_perc_{i}'] = allele_data[i-1]['rc_perc']
                sorted_row[f'Sequence_{i}'] = allele_data[i-1]['Sequence']
                sorted_row[f'UMIcount_{i}'] = allele_data[i-1]['UMIcount']
                sorted_row[f'HDist_HDRBC_{i}'] = allele_data[i-1]['HDist_HDRBC']
                sorted_row[f'AlnRefSeq_{i}'] = allele_data[i-1]['AlnRefSeq']
            else:
                sorted_row[f'Allele{i}'] = None
                sorted_row[f'rc_{i}'] = None
                sorted_row[f'rc_perc_{i}'] = None
                sorted_row[f'Sequence_{i}'] = None
                sorted_row[f'UMIcount_{i}'] = None
                sorted_row[f'HDist_HDRBC_{i}'] = None
                sorted_row[f'AlnRefSeq_{i}'] = None
        return pd.Series(sorted_row)

    # Apply the sorting function row-wise
    sorted_df = df.apply(sort_row, axis=1)

    # Merge the sorted data with the original data
    df.update(sorted_df)

    return df


# Apply the vectorized function to the DataFrame
df_sorted = sort_alleles_by_umi(df)
df = df_sorted


#classify deletions
def find_first_gap_in_alnseq(aligned_seq):
    """
    Finds the first gap position in the aligned sequence.
    """
    for i, char in enumerate(aligned_seq):
        if char == '-':
            return i
    return None

def extract_deletion_surrounding_bases(aligned_seq, gap_pos, mod_length):
    """
    Extracts bases around the gap in the aligned sequence. 
    The window extends up to the length of the deletion on both sides.
    """
    # Define the window around the gap
    start_pos = max(0, gap_pos - 25)
    
    # Extract left and right surrounding bases
    left_bases = aligned_seq[start_pos:gap_pos].replace('-', '')
    gap_end_pos = gap_pos + mod_length
    right_bases = aligned_seq[gap_end_pos:gap_end_pos + 25].replace('-', '')
    
    return left_bases, right_bases

def extract_deleted_bases(aligned_refseq, aligned_seq, mod_length):
    """
    Extracts the deleted bases from the aligned reference sequence starting from the first gap.
    """
    # Find the first gap position in the aligned sequence
    gap_pos = find_first_gap_in_alnseq(aligned_seq)
    
    if gap_pos is None:
        return None, None
    
    # Extract the reference bases where the gap starts in the aligned sequence
    deleted_bases = aligned_refseq[gap_pos:gap_pos + mod_length]
    
    return deleted_bases

def find_longest_microhomology(deleted_seq, left_part, right_part, min_homology_length=2, tolerance=2):
    """
    Finds the longest microhomology between the deleted sequence and either the left or right part, 
    using both start-to-start and end-to-end matching.
    
    Parameters:
    deleted_seq (str): The deleted sequence.
    left_part (str): The left amplicon part.
    right_part (str): The right amplicon part.
    min_homology_length (int): Minimum length of microhomology to be reported.
    tolerance (int): The number of positions within which the microhomology can start/end.
    
    Returns:
    tuple: A tuple containing the longest microhomology, its type (start-to-start or end-to-end), 
           and the source (left or right part).
    """
    
    def find_longest_start_to_start_microhomology_with_tolerance(seq1, seq2, min_homology_length=2, tolerance=2):
        longest_microhomology = ""
        longest_length = 0
        len1 = len(seq1)
        len2 = len(seq2)

        for length in range(min_homology_length, min(len1, len2) + 1):
            for start1 in range(min(tolerance + 1, len1 - length + 1)):
                for start2 in range(min(tolerance + 1, len2 - length + 1)):
                    if seq1[start1:start1 + length] == seq2[start2:start2 + length]:
                        if length > longest_length:
                            longest_microhomology = seq1[start1:start1 + length]
                            longest_length = length
        return longest_microhomology if longest_length > 0 else None

    def find_longest_end_to_end_microhomology_with_tolerance(seq1, seq2, min_homology_length=2, tolerance=2):
        longest_microhomology = ""
        longest_length = 0
        len1 = len(seq1)
        len2 = len(seq2)

        for length in range(min_homology_length, min(len1, len2) + 1):
            for end1 in range(len1 - tolerance, len1 - length + 1):
                for end2 in range(len2 - tolerance, len2 - length + 1):
                    if seq1[end1:end1 + length] == seq2[end2:end2 + length]:
                        if length > longest_length:
                            longest_microhomology = seq1[end1:end1 + length]
                            longest_length = length
        return longest_microhomology if longest_length > 0 else None

    # Find the longest start-to-start microhomology
    start_to_start_left = find_longest_start_to_start_microhomology_with_tolerance(deleted_seq, left_part, min_homology_length, tolerance)
    start_to_start_right = find_longest_start_to_start_microhomology_with_tolerance(deleted_seq, right_part, min_homology_length, tolerance)
    
    # Find the longest end-to-end microhomology
    end_to_end_left = find_longest_end_to_end_microhomology_with_tolerance(deleted_seq, left_part, min_homology_length, tolerance)
    end_to_end_right = find_longest_end_to_end_microhomology_with_tolerance(deleted_seq, right_part, min_homology_length, tolerance)
    
    # Compare the results and return the longest
    longest_microhomology = None
    microhomology_type = None
    source = None
    
    candidates = [
        (start_to_start_left, "start-to-start", "left"),
        (start_to_start_right, "start-to-start", "right"),
        (end_to_end_left, "end-to-end", "left"),
        (end_to_end_right, "end-to-end", "right")
    ]
    
    for candidate, mh_type, mh_source in candidates:
        if candidate and (not longest_microhomology or len(candidate) > len(longest_microhomology)):
            longest_microhomology = candidate
            microhomology_type = mh_type
            source = mh_source
    
    return (longest_microhomology, microhomology_type, source) if longest_microhomology else None



#classify insertions
AAV_ITR = "AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAG"
AAV_ITR_rev = "CTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTCGCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT"
GC_set = {'G', 'C'}

def has_near_match(inserted_bases, reference_sequence, max_mismatches=0):
    # Check all substrings of the reference_sequence that have the same length as inserted_bases
    for i in range(len(reference_sequence) - len(inserted_bases) + 1):
        ref_substring = reference_sequence[i:i+len(inserted_bases)]
        if hamming_distance(inserted_bases, ref_substring) <= max_mismatches:
            return True
    return False

def find_gap_in_refseq(aligned_refseq):
    # Find the position of the gap '-' in the aligned reference sequence
    gap_positions = [i for i, char in enumerate(aligned_refseq) if char == '-']
    
    # Return the first gap position found, assuming one insertion point
    if gap_positions:
        return gap_positions[0]
    else:
        return None


def extract_surrounding_bases(aligned_refseq, aligned_seq, mod_length):
    # Find the gap position in the aligned_refseq
    gap_pos = find_gap_in_refseq(aligned_refseq)
    aligned_seq = aligned_seq.rstrip('-')
    
    if gap_pos is None:
        return None, None, None
    
    # Define the start and end positions of the window around the gap (mod_length in either direction)
    start_pos = max(0, gap_pos - mod_length)
    end_pos = min(len(aligned_refseq), gap_pos + mod_length + mod_length)

    # Extract the surrounding bases around the gap from the reference and the aligned sequences
    surrounding_ref_bases = aligned_refseq[start_pos:end_pos].replace('-', '')    
    # Extract inserted bases from the aligned sequence (should be at the gap position)
    inserted_bases = str(aligned_seq[gap_pos:gap_pos + mod_length])  # Adjusted to capture `mod_length` bases
    
    return surrounding_ref_bases, inserted_bases

def find_repeating_base_element(sequence):
    # Iterate over possible substring lengths (from 1 to half the length of the sequence)
    for i in range(1, len(sequence) // 2 + 1):
        # Check if the sequence length is divisible by the current length (i)
        if len(sequence) % i == 0:
            # Get the base element (substring of length i)
            base_element = sequence[:i]
            # Repeat the base element to see if it matches the original sequence
            if base_element * (len(sequence) // i) == sequence:
                return base_element
    # If no repeating base element is found, return the original sequence
    return sequence

def classify_insertion(mod, aligned_refseq, aligned_seq):
    # Convert the modification length to an integer
    mod_length = int(mod[1:])
    
    # Extract surrounding bases and inserted bases
    surrounding_ref_bases, inserted_bases = extract_surrounding_bases(aligned_refseq, aligned_seq, mod_length)
    
    if surrounding_ref_bases is None:
        #No gap found in the aligned reference sequence
        return False
    
    # Compare the inserted bases to the surrounding sequence bases
    if inserted_bases in surrounding_ref_bases:
        #The inserted bases are within the surrounding sequence window
        return "TemplatedInsertion"
    elif len(inserted_bases) > 1 and all(base in GC_set for base in inserted_bases):
        #The inserted bases consist only of G and/or C and are longer than 1 (GC insertion).
        return "GCInsertion"
    elif len(inserted_bases) >= 8:
        if 8 <= len(inserted_bases) <= 20:
            max_mismatches = 0
        elif len(inserted_bases) > 20:
            max_mismatches = 3
        elif len(inserted_bases) > 40:
            max_mismatches = 4
        if has_near_match(inserted_bases, AAV_ITR, max_mismatches) or has_near_match(inserted_bases, AAV_ITR_rev, max_mismatches):
        # The inserted bases are AAV ITR with the allowed number of mismatches
            return "ITRInsertion"
    else:
        base_element = find_repeating_base_element(inserted_bases)
        if base_element in surrounding_ref_bases:
            return "RepeatedTemplatedInsertion"
        else:
            #The inserted bases are NOT within the surrounding sequence window or the AAV ITR.
            return False

# Function to classify the repair outcome for each allele
def classify_repair_outcomes(row, max_alleles=10):
    repair_outcomes = {}

    for i in range(1, max_alleles + 1):
        allele_column = f'Allele{i}'
        sequence_column = f'Sequence_{i}'
        alnrefseq_column = f'AlnRefSeq_{i}'
        repair_outcome_column = f'RepairOutcome_{i}'
        hamming_distance_column = f'HDist_HDRBC_{i}'
        cleaned_sequence = str(row[sequence_column]).replace('-', '')
        
        # Check if the allele exists in the row and if it has a value
        if pd.notna(row.get(allele_column)) and row[allele_column]:
            allele = row[allele_column]
            hamming_distance = row[hamming_distance_column]
            if len(allele) == hdrbc_len:
                if hamming_distance >= (hdrbc_len - 1):
                    repair_outcomes[repair_outcome_column] = 'HDR'
                    continue
                else:
                   repair_outcomes[repair_outcome_column] = f'IncompleteHDR-{int(hamming_distance)}'
                   continue
            elif allele == 'WT':
                repair_outcomes[repair_outcome_column] = "WT"
                continue
            else:
                match = re.search(r'(D(\d+))?(I(\d+))?(S(\d+))?', allele)
                if match:
                    deletion_length = int(match.group(2)) if match.group(2) else 0
                    insertion_length = int(match.group(4)) if match.group(4) else 0
                    substitution_length = int(match.group(6)) if match.group(6) else 0
                    if insertion_length > 0 and deletion_length > 0 and substitution_length > 0:
                        if HDR_anchor in cleaned_sequence:
                            repair_outcomes[repair_outcome_column] = 'IncorrectHDR'
                            continue
                        else:
                            repair_outcomes[repair_outcome_column] = "Deletion+Insertion+Substitution"
                            continue
                    elif insertion_length > 0 and deletion_length == 0 and substitution_length > 0:
                        #insertion+substitution
                        if HDR_anchor in cleaned_sequence:
                            repair_outcomes[repair_outcome_column] = 'IncorrectHDR'
                            continue
                        else:
                            mod_i = "I" + str(insertion_length)
                            insertion_class = classify_insertion(mod_i, row[alnrefseq_column], row[sequence_column])
                            if insertion_class:
                                repair_outcomes[repair_outcome_column] = str(insertion_class+"+"+"Substitution")
                            else:
                                repair_outcomes[repair_outcome_column] = "UnclassifiedInsertion+Substitution"                            
                                continue
                    elif insertion_length > 0 and deletion_length > 0 and substitution_length == 0:
                        repair_outcomes[repair_outcome_column] = "Deletion+Insertion"
                        continue
                    elif insertion_length == 0 and deletion_length > 0 and substitution_length == 0:
                        #deletions
                        deleted_seq = extract_deleted_bases(row[alnrefseq_column], row[sequence_column], deletion_length)                   
                        if pd.notna(deleted_seq):
                            left_part, right_part = extract_deletion_surrounding_bases(row[sequence_column], find_first_gap_in_alnseq(row[sequence_column]), deletion_length)
                            
                            result = find_longest_microhomology(deleted_seq, left_part, right_part)
                            if result:
                                microhomology, mh_type, source = result
                                mh_len = len(microhomology)
                                if mh_len >= 2:
                                    repair_outcomes[repair_outcome_column] = "MMEJ"
                                else:
                                    repair_outcomes[repair_outcome_column] = "NHEJ"
                            else:
                                repair_outcomes[repair_outcome_column] = "NHEJ"
                        else:
                            print(f"No Deleted seuqence detected: {row[alnrefseq_column]}")
                    elif insertion_length == 0 and deletion_length > 0 and substitution_length > 0:
                        #deletion+substitution
                        if HDR_anchor in cleaned_sequence:
                            repair_outcomes[repair_outcome_column] = 'IncorrectHDR'
                            continue
                        else:
                            deleted_seq = extract_deleted_bases(row[alnrefseq_column], row[sequence_column], deletion_length)                   
                            if pd.notna(deleted_seq):
                                left_part, right_part = extract_deletion_surrounding_bases(row[sequence_column], find_first_gap_in_alnseq(row[sequence_column]), deletion_length)
                                
                                result = find_longest_microhomology(deleted_seq, left_part, right_part)
                                if result:
                                    microhomology, mh_type, source = result
                                    mh_len = len(microhomology)
                                    
                                    if mh_len >= 2:
                                        repair_outcomes[repair_outcome_column] = "MMEJ+Substitution"
                                    else:
                                        repair_outcomes[repair_outcome_column] = "NHEJ+Substitution"
                                else:
                                    repair_outcomes[repair_outcome_column] = "NHEJ+Substitution"
                            else:
                                repair_outcomes[repair_outcome_column] = "UnclassifiedDeletion+Substitution"
                                continue
                    elif insertion_length == 0 and deletion_length == 0 and substitution_length > 0:
                        #substitution
                        if HDR_anchor in cleaned_sequence:
                            repair_outcomes[repair_outcome_column] = 'IncorrectHDR'
                            continue
                        else:
                            repair_outcomes[repair_outcome_column] = f"Substitution_{str(substitution_length)}"
                            continue
                    elif insertion_length > 0 and deletion_length == 0 and substitution_length == 0:
                        #insertion
                        insertion_class = classify_insertion(allele, row[alnrefseq_column], row[sequence_column])
                        if insertion_class:
                            repair_outcomes[repair_outcome_column] = insertion_class
                        else:
                            repair_outcomes[repair_outcome_column] = "UnclassifiedInsertion"
                        continue
                else:
                    repair_outcomes[repair_outcome_column] = ""
                    continue
        else:
            # If no value in allele column, skip processing
            repair_outcomes[repair_outcome_column] = ""
            continue
    
    return repair_outcomes

# Apply the function to each row of the dataframe
max_alleles = max(int(col.replace('Allele', '')) for col in df.columns if col.startswith('Allele'))

# Create a copy to avoid SettingWithCopyWarning
df = df.copy()

# Create a dictionary with the repair outcomes for each row and expand it into new columns
repair_outcomes_df = df.apply(lambda row: classify_repair_outcomes(row, max_alleles), axis=1, result_type='expand')

# Now concatenate the new repair outcome columns to the dataframe
df = pd.concat([df, repair_outcomes_df], axis=1)

# Rearrange the dataframe to place 'RepairOutcome_{i}' columns after 'AlnRefSeq_{i}' columns
for i in range(1, max_alleles + 1):
    alnrefseq_column = f'AlnRefSeq_{i}'
    repair_outcome_column = f'RepairOutcome_{i}'
    
    if alnrefseq_column in df.columns and repair_outcome_column in df.columns:
        # Move the repair outcome column to be immediately after the AlnRefSeq_{i} column
        cols = df.columns.tolist()
        alnrefseq_idx = cols.index(alnrefseq_column)
        cols.insert(alnrefseq_idx + 1, cols.pop(cols.index(repair_outcome_column)))
        df = df[cols]


df_sorted = df

# Create detailed repair outcome dynamically
df_sorted['DetailedRepairOutcome'] = df_sorted.apply(
    lambda row: '_'.join([row[f'Allele{i}'] for i in range(1, row['allele_count'] + 1)]),
    axis=1
)

# Function to transform HDR values
def transform_outcome(value):
    parts = value.split('_')
    new_parts = []
    for part in parts:
        if len(part) == hdrbc_len:
            new_parts.append('HDR')
        else:
            new_parts.append(part)
    return '_'.join(new_parts)

df_sorted['DetailedRepairOutcome_HDR'] = df_sorted['DetailedRepairOutcome'].apply(transform_outcome)

# Define the function to categorize alleles dynamically
def categorize_alleles(df, max_alleles):
    outcome_columns = {}
    
    for i in range(1, max_alleles + 1):
        def categorize_row(row):
            outcome_parts = row['DetailedRepairOutcome_HDR'].split('_')
            if i <= len(outcome_parts):
                return categorize_outcome(row[f'Sequence_{i}'], outcome_parts[i-1], row[f'HDist_HDRBC_{i}'])
            else:
                return None  # Return None if the allele doesn't exist
        
        outcome_columns[f'Allele{i}_Outcome'] = df.apply(categorize_row, axis=1)
    
    return pd.DataFrame(outcome_columns)

# Function to categorize repair outcome dynamically
def categorize_outcome(sequence, outcome, hamming_distance):
    category = []
    cleaned_sequence = str(sequence).replace('-', '')
    
    if pd.isna(sequence) or 'NA' in outcome or not outcome:
        return 'NotAssigned'
    if outcome == 'HDR':
        if hamming_distance >= (hdrbc_len - 1):
            return 'HDR'
        else:
           return f'IncompleteHDR-{hamming_distance}'
    if 'D' in outcome and 'S' not in outcome:
        if left_anchor not in sequence:
            category.append('PAMproximalDeletion')
        if right_anchor not in sequence:
            category.append('PAMdistalDeletion')
        if 'PAMproximalDeletion' in category and 'PAMdistalDeletion' in category:
            return 'BidirectionalDeletion'
        if HDR_anchor in cleaned_sequence:
            return 'IncorrectHDR'
    if 'I' in outcome and 'S' not in outcome:
        category.append('Insertion')
        if HDR_anchor in cleaned_sequence:
            return 'IncorrectHDR'
    if 'I' in outcome and 'D' in outcome and 'S' not in outcome:
        if HDR_anchor in cleaned_sequence:
            return 'IncorrectHDR'
    if 'S' in outcome:
        category.append('IncorrectHDR')
    
    return ', '.join(category) if category else 'WT'
    

# Apply categorization
max_alleles = df_sorted['allele_count'].max()
df_outcomes = categorize_alleles(df_sorted, max_alleles)

# Combine outcomes into a new column
df_sorted['Combined_Outcome'] = df_outcomes.apply(
    lambda row: '_'.join(filter(None, [row[f'Allele{i}_Outcome'] for i in range(1, max_alleles + 1)])),
    axis=1
)

# Create HDR columns and concatenate them
hdr_columns = {
    f'Allele{i}_HDR': df_sorted[f'Allele{i}'].apply(lambda x: 'HDR' if len(str(x)) == hdrbc_len else x)
    for i in range(1, max_alleles + 1)
}
df_hdr = pd.DataFrame(hdr_columns)

# Concatenate all new columns at once to the original DataFrame
df_combined = pd.concat([df_sorted, df_outcomes, df_hdr], axis=1)

# De-fragment the DataFrame by creating a copy
df_combined = df_combined.copy()

# Frequency analysis
allele_freq = pd.Series(dtype=float)
for i in range(1, df_combined['allele_count'].max() + 1):
    allele_freq = allele_freq.add(df_combined[f'RepairOutcome_{i}'].value_counts(), fill_value=0)

pd.DataFrame(
    {"Allele": allele_freq.index, "Frequency": allele_freq.values}
).to_excel("AlleleFrequencies.xlsx", index=False)
plt.figure(figsize=(40, 6))
plt.bar(allele_freq.index, allele_freq.values, color='purple')
plt.title('Combined Frequency of All Alleles')
plt.xlabel('Allele')
plt.ylabel('Frequency')
plt.savefig('allele_frequencies.png', format='png', dpi=300)
plt.close()

# Categorization of deletions and insertions dynamically
def categorize_deletions(deletion):
    try:
        size = int(deletion[1:])
        if size <= 8:
            return 'D1-D8'
        elif 9 <= size <= 30:
            return 'D9-D30'
        else:
            return 'D31+'
    except ValueError:
        if deletion == 'HDR':
            return 'HDR'
        elif deletion == 'WT':
            return 'WT'
        else:
            return 'NA'
    except TypeError:
        return 'NA'

def categorize_insertions(insertion):
    try:
        size = int(insertion[1:])
        if size <= 3:
            return 'I1-I3'
        elif 4 <= size <= 20:
            return 'I4-I20'
        else:
            return 'I21+'
    except ValueError:
        if insertion == 'HDR':
            return 'HDR'
        elif insertion == 'WT':
            return 'WT'
        else:
            return 'NA'
    except TypeError:
        return 'NA'

# Apply categorization dynamically
sized_columns = {
    f'Allele{i}_Outcome_sized': df_combined[f'Allele{i}_HDR'].apply(
        lambda x: categorize_deletions(x) if 'D' in str(x) else (categorize_insertions(x) if 'I' in str(x) else x)
    )
    for i in range(1, df_combined['allele_count'].max() + 1)
}

# Concatenate new sized outcome columns to the DataFrame
df_sized = pd.concat([df_combined, pd.DataFrame(sized_columns)], axis=1)

# De-fragment the DataFrame by creating a copy
df_sized = df_sized.copy()

# Frequency analysis for sized outcomes
allele_freq_sized = pd.Series(dtype=float)
for i in range(1, df_sized['allele_count'].max() + 1):
    allele_freq_sized = allele_freq_sized.add(df_sized[f'Allele{i}_Outcome_sized'].value_counts(), fill_value=0)

# Save the frequency analysis and the updated DataFrame to Excel files
pd.DataFrame(
    {"Allele": allele_freq_sized.index, "Frequency": allele_freq_sized.values}
).to_excel("AlleleFrequencies_sized.xlsx", index=False)
df_sized.to_excel('EditingOutcomesByCellBarcode.xlsx', index=False, engine='openpyxl')

# Determine the number of rows and columns for the subplots grid based on the maximum number of alleles
n_cols = 2  # Keep 2 columns
n_rows = math.ceil(max_alleles / n_cols)  # Calculate the number of rows needed

# Plotting outcome counts for each allele
plt.figure(figsize=(200, 150), constrained_layout=True)
for i in range(1, max_alleles + 1):
    plt.subplot(n_rows, n_cols, i)
    outcome_counts = df_sized[f'Allele{i}_Outcome'].value_counts()
    if not outcome_counts.empty:
        outcome_counts.plot(kind='bar')
        plt.title(f'Allele{i} Outcome Counts')
        plt.xlabel('Outcome Categories')
        plt.ylabel('Counts')

plt.savefig('allele_outcome_counts.png', format='png', dpi=300)
plt.close()

# Plotting combined outcome counts
plt.figure(figsize=(200, 150), constrained_layout=True)
df_sized['Combined_Outcome'].value_counts().plot(kind='bar')
plt.title('Combined Outcome Counts')
plt.xlabel('Outcome Categories')
plt.ylabel('Counts')
plt.savefig('allele_combinedoutcome_counts.png', format='png', dpi=300)
plt.close()

cleanup_dro_outputs()
