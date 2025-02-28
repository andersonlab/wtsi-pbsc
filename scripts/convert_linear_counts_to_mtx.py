import pandas as pd
import argparse


def create_mtx(isoquant_counts_dat):


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

    return [full_mtx,features_dat.drop(columns=['feature_idx']),groups_dat.drop(columns=['group_idx'])]










def write_mtx_dat(full_mtx,features_dat,groups_dat,output_dir,output_prefix):
        if output_prefix.strip() != '':
            output_prefix=output_prefix+'.'
        #First line data
        firstline_dat=' '.join([str(features_dat.shape[0]), str(groups_dat.shape[0]),str(full_mtx.shape[0])])+'\n'

        #Headers and mtx file
        full_mtx_f=open(output_dir+output_prefix+'matrix.mtx','w')
        full_mtx_f.write("%%MatrixMarket matrix coordinate real general\n")
        full_mtx_f.write("%\n")
        full_mtx_f.write(firstline_dat)
        full_mtx_f.close()
        full_mtx_f=open(output_dir+output_prefix+'matrix.mtx','a')
        full_mtx.to_csv(full_mtx_f,sep=' ',header=None,index=False)
        full_mtx_f.close()
        print('wrote matrix.mtx file...')

        #genes.tsv file
        features_dat[1]=features_dat['#feature_id']
        features_dat.to_csv(output_dir+output_prefix+'genes.tsv',sep='\t',header=False,index=False)
        print('wrote genes.tsv file...')

        #barcodes.tsv file
        groups_dat.to_csv(output_dir+output_prefix+'barcodes.tsv',sep='\t',header=False,index=False)
        print('wrote barcodes.tsv file...')
        print('==================')

def load_all_isoquant_dat(isoquant_count_fs):
    all_isoqunat_counts_dat=pd.DataFrame()
    for isoquant_count_f in isoquant_count_fs:
        isoquant_counts_dat=pd.read_csv(isoquant_count_f,sep='\t')
        isoquant_counts_dat=isoquant_counts_dat.loc[isoquant_counts_dat['count']!=0,:]
        all_isoqunat_counts_dat=pd.concat([all_isoqunat_counts_dat,isoquant_counts_dat],ignore_index=True)
        print('Loaded isoquant (linear) counts from {f}...'.format(f=isoquant_count_f))
    return all_isoqunat_counts_dat


def main():
    parser = argparse.ArgumentParser(description="Process Isoquant count file with specified parameters.")

    parser.add_argument("-i",'--isoquant_files', dest='isoquant_count_fs',nargs='+',required=True,type=str, help="Path to the Isoquant count file.")
    parser.add_argument('-d',"--output_dir", dest='output_dir',default='',type=str, help="Directory where output files will be saved.")
    parser.add_argument('-p',"--output_prefix", dest='output_prefix',default='',type=str, help="Prefix for output file names.")

    args = parser.parse_args()


    all_mtx_dat=pd.DataFrame()
    all_features_dat=pd.DataFrame()
    all_groups_dat=pd.DataFrame()

    output_prefix=args.output_prefix
    output_dir=args.output_dir
    isoquant_count_fs=args.isoquant_count_fs



    all_isoqunat_counts_dat=load_all_isoquant_dat(isoquant_count_fs)
    print('=========')
    print('Loaded all isoquant files')
    print('=========')
    full_mtx,features_dat,groups_dat=create_mtx(all_isoqunat_counts_dat)
    print('=========')
    print('Created MTX files')
    print('=========')
    write_mtx_dat(full_mtx,features_dat,groups_dat,output_dir,output_prefix)
    print('=========')
    print('Wrote MTX files')
    print('=========')
if __name__ == "__main__":
    main()
