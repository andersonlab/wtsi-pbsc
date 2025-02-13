import pandas as pd
import argparse


def create_mtx(isoquant_f):
    isoquant_counts_dat=pd.read_csv(isoquant_f,sep='\t')
    print('Loaded isoquant (linear) counts...')

    #Extracting feature index
    features_dat=isoquant_counts_dat['#feature_id'].drop_duplicates().reset_index(drop=True).reset_index().rename(columns={'index':'feature_idx'})
    features_dat['feature_idx']=features_dat['feature_idx']+1
    print('Generated genes index...')


    #Extracting cell index
    groups_dat=isoquant_counts_dat['group_id'].drop_duplicates().reset_index(drop=True).reset_index().rename(columns={'index':'group_idx'})
    groups_dat['group_idx']=groups_dat['group_idx']+1
    print('Generated barcodes index...')


    #Creating Full mtx dataframe
    full_mtx=(isoquant_counts_dat.merge(features_dat,on='#feature_id',how='left').merge(groups_dat,on='group_id',how='left')).loc[:,['feature_idx','group_idx','count']]
    #Removing 0s to create sparse
    full_mtx=full_mtx.loc[full_mtx['count']!=0,:]
    print('Filtered out 0 counts...')

    return [features_dat,groups_dat,full_mtx]










def write_mtx_dat(full_mtx,features_dat,groups_dat,output_dir,output_prefix):
        #First line data
        firstline_dat=' '.join([str(features_dat.shape[0]), str(groups_dat.shape[0]),str(full_mtx.shape[0])])+'\n'

        #Headers and mtx file
        full_mtx_f=open(output_dir+'matrix.mtx','w')
        full_mtx_f.write("%%MatrixMarket matrix coordinate real general\n")
        full_mtx_f.write("%\n")
        full_mtx_f.write(firstline_dat)
        full_mtx_f.close()
        full_mtx_f=open(output_dir+output_prefix+'matrix.mtx','a')
        full_mtx.to_csv(full_mtx_f,sep=' ',header=None,index=False)
        full_mtx_f.close()
        print('wrote matrix.mtx file...')

        #genes.tsv file
        features_dat.to_csv(output_dir+output_prefix+'genes.tsv',sep='\t',header=False,index=False)
        print('wrote genes.tsv file...')

        #barcodes.tsv file
        groups_dat.to_csv(output_dir+output_prefix+'barcodes.tsv',sep='\t',header=False,index=False)
        print('wrote barcodes.tsv file...')
        print('==================')


def main():
    parser = argparse.ArgumentParser(description="Process Isoquant count file with specified parameters.")

    parser.add_argument("-i",'--isoquant_files', dest='isoquant_count_fs',nargs='+',required=True,type=str, help="Path to the Isoquant count file.")
    parser.add_argument('-d',"--output_dir", dest='output_dir',required=True,type=str, help="Directory where output files will be saved.")
    parser.add_argument('-p',"--output_prefix", dest='output_prefix',required=True,type=str, help="Prefix for output file names.")

    args = parser.parse_args()


    all_mtx_dat=pd.DataFrame()
    all_features_dat=pd.DataFrame()
    all_groups_dat=pd.DataFrame()

    for isoquant_f in args.isoquant_fs
        full_mtx,features_dat,groups_dat=create_mtx(isoquant_f)
        all_mtx_dat=pd.concat([all_mtx_dat,full_mtx],ignore_index=True)
        all_features_dat=pd.concat([all_features_dat,features_dat],ignore_index=True)
        all_groups_dat=pd.concat([all_groups_dat,groups_dat],ignore_index=True)

    write_mtx_dat(all_mtx_dat,all_features_dat,all_groups_dat,args.output_dir,args.output_prefix)
if __name__ == "__main__":
    main()
