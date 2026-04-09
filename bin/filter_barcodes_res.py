#!/usr/bin/env python3
import pandas as pd
import sys

def filter_with_merge(input_file, threshold1, threshold2, output_file):
    # Read data
    df = pd.read_csv(input_file, sep='\t')
    
    # Step 1: Get unique combinations (preserves original order)
    all_combinations = df[['pool', 'donor']].drop_duplicates()
    print(f"Unique (pool, donor) combinations: {len(all_combinations)}")
    
    # Step 2: Filter by both thresholds
    filtered_df = df[
        (df['n_barcodes_in_common'] >= int(threshold1)) & 
        (df['%_in_common'] >= float(threshold2))
    ]
    print(f"Combinations passing filters (n_barcodes_in_common>={threshold1} & %_in_common>={threshold2}): {len(filtered_df)}")
    
    # Step 3: Check for multiple matches per (pool, donor) combination
    # Group by pool and donor, find max %_in_common
    grouped = filtered_df.groupby(['pool', 'donor'])['%_in_common'].max().reset_index()
    grouped.columns = ['pool', 'donor', 'max_pct']
    
    # Merge to identify rows with max %_in_common
    filtered_with_max = filtered_df.merge(grouped, on=['pool', 'donor'])
    filtered_with_max = filtered_with_max[filtered_with_max['%_in_common'] == filtered_with_max['max_pct']]
    
    # Step 4: Identify combinations with multiple rows having the same max %_in_common
    count_per_combination = filtered_with_max.groupby(['pool', 'donor']).size().reset_index(name='count')
    multiple_matches = count_per_combination[count_per_combination['count'] > 1]
    
    if len(multiple_matches) > 0:
        print(f"⚠ Found {len(multiple_matches)} combinations with multiple matches at highest %_in_common")
        
        # Get the actual rows with multiple matches
        multiple_rows = filtered_with_max.merge(
            multiple_matches[['pool', 'donor']], 
            on=['pool', 'donor']
        )
        
        # Save to multiple_good_matches.tsv
        multiple_rows_output = multiple_rows.drop(columns=['max_pct'])
        multiple_rows_output.to_csv('multiple_good_barcode_matches.tsv', sep='\t', index=False)
        print(f"✓ Multiple matches saved to: multiple_good_barcode_matches.tsv ({len(multiple_rows)} rows)")
        
        # Remove these combinations from filtered_with_max
        filtered_with_max = filtered_with_max.merge(
            multiple_matches[['pool', 'donor']], 
            on=['pool', 'donor'], 
            how='left', 
            indicator=True
        )
        filtered_with_max = filtered_with_max[filtered_with_max['_merge'] == 'left_only']
        filtered_with_max = filtered_with_max.drop(columns=['_merge'])
    
    # Step 5: Keep only necessary columns (drop max_pct helper column)
    filtered_with_max = filtered_with_max.drop(columns=['max_pct'])
    
    # Step 6: Left merge with all combinations
    result = all_combinations.merge(
        filtered_with_max,
        on=['pool', 'donor'],
        how='left'
    )
    
    # Step 7: Fill missing values
    fill_values = {
        'barcode_sample_id': 'NONE',
        'n_barcodes_in_common': 0,
        'n_barcodes': 0,
        '%_in_common': 0.0,
        'rank': 0
    }
    result = result.fillna(fill_values)
    
    # Convert to appropriate types
    result['n_barcodes_in_common'] = result['n_barcodes_in_common'].astype(int)
    result['n_barcodes'] = result['n_barcodes'].astype(int)
    result['rank'] = result['rank'].astype(int)
    
    # Save
    result.to_csv(output_file, index=False)
    
    print(f"Combinations filled with NONE/0: {len(all_combinations) - len(filtered_with_max)}")
    print(f"\n✓ Output saved to: {output_file} ({len(result)} rows)")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: filter_barcodes_res.py input.tsv threshold1 threshold2 output.csv")
        print("  threshold1: minimum number of common barcodes (n_barcodes_in_common)")
        print("  threshold2: minimum percentage in common (%_in_common)")
        sys.exit(1)
    
    filter_with_merge(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])