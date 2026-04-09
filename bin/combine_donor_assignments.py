#!/usr/bin/env python3
import pandas as pd
import sys
import os

def process_donor_files(file_list, table1_file, table2_file, output_file):
    
    # Step 1: Read the file list
    print(f"Reading file list: {file_list}")
    file_mapping = pd.read_csv(file_list, sep='\t', header=None, names=['pool', 'file_path'])
    print(f"Found {len(file_mapping)} pools")
    
    # Step 2: Extract unique donors from each donor_ids.tsv
    all_donors = []
    
    for idx, row in file_mapping.iterrows():
        pool_id = row['pool']
        file_path = row['file_path']
        
        try:
            # Read donor_ids.tsv file
            donor_df = pd.read_csv(file_path, sep='\t')
            
            # Get unique donor_ids
            unique_donors = donor_df['donor_id'].unique()
            
            # Create records for this pool
            for donor in unique_donors:
                all_donors.append({'pool': pool_id, 'donor': donor})
            
            print(f"  {pool_id}: {len(unique_donors)} unique donors")
            
        except Exception as e:
            print(f"  WARNING: Could not read {file_path}: {e}")
            continue
    
    # Create base dataframe
    base_df = pd.DataFrame(all_donors)
    base_df=base_df[(base_df['donor'] != 'doublet') & (base_df['donor'] != 'unassigned')]
    
    # Step 3: Read table1 and merge donor_gt
    if os.path.exists(table1_file):
        print(f"\nReading table1: {table1_file}")
        table1 = pd.read_csv(table1_file).rename(columns={'donor_query': 'donor'})

        base_df = base_df.merge(
            table1[['pool', 'donor', 'donor_gt']],
            on=['pool', 'donor'],
            how='left'
        )
        print(f"  Merged donor_gt: {base_df['donor_gt'].notna().sum()} matches")
    
    # Step 4: Read table2 and merge barcode_sample_id
    if os.path.exists(table2_file):
        print(f"\nReading table2: {table2_file}")
        table2 = pd.read_csv(table2_file)

        base_df = base_df.merge(
            table2[['pool', 'donor', 'barcode_sample_id']],
            on=['pool', 'donor'],
            how='left'
        )
        print(f"  Merged barcode_sample_id: {base_df['barcode_sample_id'].notna().sum()} matches")
    
    # Step 5: Save result
    base_df.to_csv(output_file, index=False)
    

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python process_donors.py file_list.tsv, gt_res.csv barcode_res.csv output.csv")
        print("\nfile_list.txt format:")
        print("  POOL_ID<tab>PATH_TO_DONOR_IDS_TSV")
        sys.exit(1)
    
    process_donor_files(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])