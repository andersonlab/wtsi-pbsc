import pandas as pd
import argparse
import os


def is_valid_counts_file(file_path):
    """Check if the file contains non-comment lines."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return False
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):  # Check for at least one non-comment line
                return True
    return False

if __name__=="__main__":





    parser = argparse.ArgumentParser(description="Process intron/exon grouped counts.")

    parser.add_argument('-i',"--input_f", dest='input_f', type=str, required=True, help="Input file (TSV format).")
    parser.add_argument('-o',"--output_dir", dest='output_dir', type=str, default="", help="Output directory.")
    parser.add_argument('-p',"--output_prefix", dest='output_prefix', type=str, default="OUT", help="Prefix for output files.")

    args = parser.parse_args()
    output_dir=args.output_dir
    output_prefix=args.output_prefix
    input_f=args.input_f

    if is_valid_counts_file(input_f):
        input_dat=pd.read_csv(input_f,sep='\t')

        input_dat['#feature_id']=input_dat['#chr'].astype(str)+":"+\
        input_dat['start'].astype(str)+":"+\
        input_dat['end'].astype(str)+":"+\
        input_dat['strand'].astype(str)+":"+\
        input_dat['flags'].astype(str)+"_"+\
        input_dat['gene_ids'].astype(str)


        input_dat_include_counts=input_dat[['#feature_id','group_id','include_counts']].rename(columns={'include_counts': 'count'})
        input_dat_exclude_counts=input_dat[['#feature_id','group_id','exclude_counts']].rename(columns={'exclude_counts': 'count'})
    else:
        input_dat_include_counts=pd.DataFrame({'#feature_id':[],'group_id':[],'count':[]})
        input_dat_exclude_counts=pd.DataFrame({'#feature_id':[],'group_id':[],'count':[]})
    input_dat_include_counts.to_csv(output_dir+output_prefix+'.include_counts.tsv',sep='\t',index=False)
    input_dat_exclude_counts.to_csv(output_dir+output_prefix+'.exclude_counts.tsv',sep='\t',index=False)
