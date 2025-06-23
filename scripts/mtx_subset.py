import pandas as pd
from convert_linear_counts_to_mtx import create_mtx
from convert_linear_counts_to_mtx import write_mtx_dat
import argparse


def mtx_subset(mtx_dat,genes_dat,barcodes_dat,subset_list):
    mtx_dat=mtx_dat.rename(columns={0:'feature_idx',1:'barcode_idx',2:'count'})
    genes_dat=genes_dat.loc[:,[0]].reset_index().rename(columns={'index':'feature_idx',0:'#feature_id'})
    genes_dat['feature_idx']=genes_dat['feature_idx']+1

    barcodes_dat=barcodes_dat.loc[:,[0]].reset_index().rename(columns={'index':'barcode_idx',0:'group_id'})
    barcodes_dat['barcode_idx']=barcodes_dat['barcode_idx']+1
    full_mtx_dat=mtx_dat.merge(genes_dat,how='left',on='feature_idx').merge(barcodes_dat,how='left',on='barcode_idx')

    n_na_features=full_mtx_dat.loc[full_mtx_dat['#feature_id'].isna(),:].shape[0]
    n_na_barcodes=full_mtx_dat.loc[full_mtx_dat['group_id'].isna(),:].shape[0]
    print("Number of features in matrix.mtx with no matching id in genes.tsv: {n}".format(n=n_na_features))
    print("Number of barcodes in matrix.mtx with no matching id in barcodes.tsv: {n}".format(n=n_na_barcodes))

    if n_na_features > 0:
        raise ValueError('There are non-matching features between MTX and genes.tsv')
    if n_na_barcodes > 0:
        raise ValueError('There are non-matching barcodes between MTX and barcodes.tsv')
    full_mtx_dat=full_mtx_dat.loc[full_mtx_dat['#feature_id'].isin(subset_list),:]
    new_mtx_dat,new_genes_dat,new_barcodes_dat=create_mtx(full_mtx_dat.drop(columns=['feature_idx','barcode_idx']) )
    return [new_mtx_dat,new_genes_dat,new_barcodes_dat]

def main():

    parser = argparse.ArgumentParser(description="Subset MTX (Needs create_mtx and write_mtx from convert_linear_counts_to_mtx.py)")

    parser.add_argument("-i",'--input_dir', dest='input_dir',required=True,type=str, help="Path to input MTX directory.")
    parser.add_argument("-s",'--subset_list', dest='subset_list_f',required=True,type=str, help="Feature subset list.")

    parser.add_argument('-d',"--output_dir", dest='output_dir',default='',type=str, help="Directory where output files will be saved.")
    parser.add_argument('-p',"--output_prefix", dest='output_prefix',default='',type=str, help="Prefix for output file names.")

    args = parser.parse_args()




    input_dir=args.input_dir
    output_dir=args.output_dir
    output_prefix=args.output_prefix
    subset_list_f=args.subset_list_f

    mtx_f=input_dir+'matrix.mtx'
    genes_f=input_dir+'genes.tsv'
    barcodes_f=input_dir+'barcodes.tsv'

    mtx_dat=pd.read_csv(mtx_f,sep=' ',skiprows=3,header=None)
    genes_dat=pd.read_csv(genes_f,sep='\t',header=None)
    barcodes_dat=pd.read_csv(barcodes_f,sep='\t',header=None)
    subset_list_dat=pd.read_csv(subset_list_f,sep='\t',header=None)
    subset_list=subset_list_dat[0].tolist()

    new_mtx_dat,new_genes_dat,new_barcodes_dat=mtx_subset(mtx_dat,genes_dat,barcodes_dat,subset_list)
    write_mtx_dat(new_mtx_dat,new_genes_dat,new_barcodes_dat,output_dir,output_prefix)


if __name__ == "__main__":
    main()
