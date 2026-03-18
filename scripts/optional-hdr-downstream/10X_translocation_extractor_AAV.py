# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 13:19:08 2024

@author: moritzschlapansky
"""

import pandas as pd
import sys
from collections import defaultdict

translocation_file = sys.argv[1]

translocation_df = pd.read_csv(translocation_file, header=0, low_memory=False)
translocation_df = translocation_df.fillna(0)

print(f"Number of unique barcodes before filtering: {translocation_df['bc'].nunique()}")
# Filter out rows where 'chromosome' is '*' (unaligned reads)
unaligned_translocations_df = translocation_df[translocation_df['chromosome'] == '*']
unaligned_translocations_df.to_csv('unaligned_translocations.csv', index=False)
translocation_df = translocation_df[translocation_df['chromosome'] != '*']
translocation_df = translocation_df[translocation_df['read_count'] >= 3]

print(f"Number of unique barcodes after filtering: {translocation_df['bc'].nunique()}")

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))




# Define a function to merge overlapping or close locations
def merge_locations(df, distance_threshold=150):
    merged = []
    current_start = df.iloc[0]['start']
    current_end = df.iloc[0]['end']
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

output_columns = ['bc', 'total_umis', 'total_reads', 'total_unique_locations']
output_file = translocation_file.rstrip('.csv') + "_output.csv"

if not output:
    # Downstream remapping can legitimately produce no passing barcodes.
    # In that case, write an empty table with the base columns instead of failing.
    pd.DataFrame(columns=output_columns).to_csv(output_file, index=False)
    print("No AAV translocations passed filtering; wrote an empty output table.")
    sys.exit(0)

# Convert to a DataFrame for better readability
max_locations = max(len(result) for result in output) - 4  # to account for the first 4 columns

# Add dynamic column headers based on the number of locations
for i in range(max_locations // 4):  # every location has 4 columns (location, umi_count, read_count, read_seq)
    output_columns.extend([f'location_{i+1}', f'umi_count_loc_{i+1}', f'read_count_loc_{i+1}', f'read_seq_loc_{i+1}'])

output_df = pd.DataFrame(output, columns=output_columns)
output_df.to_csv(output_file, index=False)
