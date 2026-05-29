import pandas as pd
import argparse

def filter_bam_stats_barcodes(barcodes_file,min_umi=1000):
    barcodes_dat=pd.read_csv(barcodes_file,sep="\t")
    barcodes_dat=barcodes_dat.loc[(barcodes_dat["NumberOfMolecules"]>=min_umi) | (barcodes_dat["RealCell"]=="cell"),:] 
    return barcodes_dat.loc[:,["#BarcodeSequence"]].drop_duplicates().reset_index(drop=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter barcodes based on UMI counts.')
    parser.add_argument('-b','--barcodes_file', type=str, help='Path to the barcodes file (tab-separated).')
    parser.add_argument('-o','--output_file', type=str, help='Path to the output file (tab-separated).')
    parser.add_argument('--min_umi', type=int, default=1000, help='Minimum number of UMIs required to keep a barcode (default: 1000).')
    
    args = parser.parse_args()
    
    filtered_barcodes = filter_bam_stats_barcodes(args.barcodes_file, args.min_umi)
    filtered_barcodes.to_csv(args.output_file, sep="\t", index=False,header=False)
