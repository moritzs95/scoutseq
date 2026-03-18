# -*- coding: utf-8 -*-
"""
Summarize remapped PARSE translocation coordinates per cell barcode.
"""

import pandas as pd
import sys
from collections import defaultdict

sub_library = ""
cbc_path = sys.argv[1]
translocation_file = sys.argv[2]

translocation_df = pd.read_csv(translocation_file, header=0, low_memory=False)
translocation_df = translocation_df.fillna(0)

print(f"Number of unique barcodes before filtering: {translocation_df['bc'].nunique()}")
# Filter out rows where 'chromosome' is '*' (unaligned reads)
unaligned_translocations_df = translocation_df[translocation_df['chromosome'] == '*']
unaligned_translocations_df.to_csv('unaligned_translocations.csv', index=False)
translocation_df = translocation_df[translocation_df['chromosome'] != '*']
translocation_df = translocation_df[translocation_df['read_count'] >= 3]



def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

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
    translocation_df = reformat_barcodes(translocation_df, 'bc', sub_library)

# Read the cell metadata
cell_metadata = pd.read_csv(f"{cbc_path}cell_metadata_unfiltered.csv", header=0, low_memory=False)
whitelist_unfiltered = set(cell_metadata['bc_wells'].tolist())
cell_metadata.set_index('bc_wells', inplace=True)


all_barcodes_before_filtering = set(translocation_df['bc'].unique())

whitelist_df = pd.read_csv(f"{cbc_path}cbcs_allsamples_filtered.csv", header=0, low_memory=False)
whitelist = set(whitelist_df['Cells.all_samples.data.'].tolist())

barcodes_filtered = set(whitelist)
filtered_overlap = all_barcodes_before_filtering.intersection(barcodes_filtered)
num_filtered_overlap = len(filtered_overlap)

# Sort rows so that rows containing '-' in 'bc' column are moved to the end
translocation_df = translocation_df.sort_values(by='bc', key=lambda col: col.str.contains('-'))

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
translocation_df_updated = update_dataframe_barcodes_optimized_trie(translocation_df, 'bc', whitelist_unfiltered)

#remove NA cells
translocation_df_updated = translocation_df_updated[translocation_df_updated['bc'].notna()]


# Get the count of unique values in these filtered rows
rows_with_hyphen = translocation_df_updated[translocation_df_updated['bc'].str.contains('-')]
unique_count = rows_with_hyphen['bc'].nunique()
print("Number of unique rows containing incomplete barcode:", unique_count)

translocation_df = translocation_df_updated



# Consolidate identical UMIs (can improve this part)
translocation_df = translocation_df.groupby(['bc', 'umi']).agg({
    'read_seq': 'first',
    'mod_translocation': 'first',
    'read_count': 'sum',
    'chromosome': 'first',
    'start': 'first',
    'end': 'first',
    'location': 'first',
}).reset_index()

print(f"Number of unique barcodes after filtering: {translocation_df['bc'].nunique()}")
# Define a function to merge overlapping or close locations
def merge_locations(df, distance_threshold=150):
    merged = []
    current_start = int(df.iloc[0]['start'])
    current_end = int(df.iloc[0]['end'])
    current_chr = df.iloc[0]['chromosome']
    
    for i in range(1, len(df)):
        row = df.iloc[i]
        if row['chromosome'] == current_chr and row['start'] <= current_end + distance_threshold:
            # If the current location is close enough, merge by updating the end
            current_end = max(current_end, row['end'])
        else:
            # Otherwise, add the current merged location to the result and move to the next
            merged.append([current_chr, current_start, current_end])
            current_start = row['start']
            current_end = row['end']
            current_chr = row['chromosome']
    
    # Don't forget to add the last merged region
    merged.append([current_chr, current_start, current_end])
    
    return pd.DataFrame(merged, columns=['chromosome', 'start', 'end'])

# Group by 'bc' and summarize the information
output = []

for bc, group in translocation_df.groupby('bc'):
    total_umis = group['umi'].nunique()
    total_reads = group['read_count'].sum()

    # Sort the group by chromosome and start position for location merging
    group = group.sort_values(by=['chromosome', 'start'])
    
    # Merge locations that are within 150 bases of each other
    merged_locations = merge_locations(group[['chromosome', 'start', 'end']])
    
    # Add information about the total number of unique locations after merging
    total_unique_locations = len(merged_locations)

    # Initialize the base structure for this bc
    result = [bc, total_umis, total_reads, total_unique_locations]
    
    # Add information about each merged location
    for i, row in merged_locations.iterrows():
        loc_group = group[(group['chromosome'] == row['chromosome']) & 
                          (group['start'] >= row['start']) & 
                          (group['end'] <= row['end'])]
        location = f"{row['chromosome']}:{row['start']}-{row['end']}"
        umi_count_loc = loc_group['umi'].nunique()
        read_count_loc = loc_group['read_count'].sum()
        read_seq_loc = loc_group['read_seq'].iloc[0]  # assuming the same read_seq for the same location

        # Append the location details
        result.extend([location, umi_count_loc, read_count_loc, read_seq_loc])
    
    # Add the result to the output list
    output.append(result)

# Convert to a DataFrame for better readability
output_columns = ['bc', 'total_umis', 'total_reads', 'total_unique_locations']
max_locations = max(len(result) for result in output) - 4  # to account for the first 4 columns

# Add dynamic column headers based on the number of locations
for i in range(max_locations // 4):  # every location has 4 columns (location, umi_count, read_count, read_seq)
    output_columns.extend([f'location_{i+1}', f'umi_count_loc_{i+1}', f'read_count_loc_{i+1}', f'read_seq_loc_{i+1}'])

# Create a new DataFrame
output_df = pd.DataFrame(output, columns=output_columns)
output_file = translocation_file.rstrip('.csv') + "_output.csv"

# Display the result
output_df.to_csv(output_file, index=False)
