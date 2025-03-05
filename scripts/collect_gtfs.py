def is_valid_gtf(file_path):
    """Check if the GTF file contains non-comment lines."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return False
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):  # Check for at least one non-comment line
                return True
    return False
def process_ref_db(ref_db_f,ref_gtf_f):
    ref_db=[]
    if ref_db_f is not None:
        ref_db=[gffutils.FeatureDB(ref_db_f,keep_order=True)]
    else:
        if ref_gtf_f is not None:
            ref_db=create_feature_db_from_gtf([ref_gtf_f])
    return(ref_db)
def create_feature_db_from_gtf(query_gtf_fs):
    query_dbs=[]
    for query_gtf_f in query_gtf_fs:
        if is_valid_gtf(query_gtf_f):
            query_db = gffutils.create_db(
                query_gtf_f, dbfn=":memory:", force=True, keep_order=True,
                disable_infer_transcripts=True, disable_infer_genes=True
            )
            query_dbs.append(query_db)
            print('DB created for {0}'.format(query_gtf_f))
        else:
            print('WARNING: Skipping empty GTF file:', query_gtf_f)


    return(query_dbs)

def clean_record_attributes(record):
    cleaned_record_attributes = {k: v for k, v in record.attributes.items() if v and v[0] != ""}
    for att in cleaned_record_attributes:

        formatted_record_atts=[ att_val.replace(';', '').replace('"', '') for att_val in cleaned_record_attributes[att] ]
        cleaned_record_attributes[att]=formatted_record_atts
    return(cleaned_record_attributes)

def recount_gene_transcripts(gene_dict):
    for gene_id in gene_dict.keys():
        gene_dict[gene_id]['gene'].attributes['transcripts']=[str(len(gene_dict[gene_id]['transcripts']))]
    return(gene_dict)

def merge_dbs(query_dbs):
    gene_dict={}
    for i,query_db in enumerate(query_dbs):
        print("Processing {x}st DB...".format(x=i+1))

        for j,gene in enumerate(query_db.features_of_type("gene")):
            gene_id = gene.attributes.get("gene_id")[0]

            if gene_id is not None:
                if gene_id not in gene_dict:
                    gene.attributes=clean_record_attributes(gene)
                    gene_dict[gene_id] = {"gene":gene,"transcripts": {}}

                # Add children (transcripts, exons, CDS, etc.)
                for transcript in query_db.children(gene,featuretype="transcript"):
                    transcript_id=transcript.attributes.get("transcript_id")[0]

                    if transcript_id in gene_dict[gene_id]['transcripts'].keys():
                        print('WARNING: duplicate transcript {0}. Used first instance'.format(transcript_id))
                        continue
                    transcript.attributes=clean_record_attributes(transcript)
                    for child in query_db.children(transcript):
                        child.attributes=clean_record_attributes(child)
                        if child.featuretype != 'transcript':
                            if transcript_id in gene_dict[gene_id]['transcripts'].keys():
                                gene_dict[gene_id]['transcripts'][transcript_id]['children'].append(child)
                            else:
                                gene_dict[gene_id]['transcripts'][transcript_id]={'transcript':transcript,'children':[child]}

    gene_dict=recount_gene_transcripts(gene_dict)

    return(gene_dict)

def format_gtf_output(gene_dict):
    output_str=""
    for gene,info in gene_dict.items():
        output_str+=str(info['gene']).rstrip(';')+";\n"
        for transcript,transcript_info in info['transcripts'].items():
            output_str+=str(transcript_info['transcript']).rstrip(';')+";\n"
            for child in transcript_info['children']:
                output_str+=str(child).rstrip(';')+";\n"
    return(output_str)

#TODO: enable subsetting of reference
def subset_db(db,transcript_subset):
    pass




import gffutils
import pandas as pd
import argparse
import os





def main():
    parser = argparse.ArgumentParser(description="Process multiple GTF files with a reference GTF and selected transcripts.")

    parser.add_argument(
        "-q", "--query_gtf_files",
        dest="query_gtf_fs",
        nargs="+",  # Allows multiple files
        required=False,
        default=None,
        help="Input query GTF files (multiple allowed)"
    )
    parser.add_argument(
        "-Q", "--query_gtf_fofn",
        dest="query_gtf_fofn",
        required=False,
        default=None,
        help="Path to file listing input query GTF files"
    )

    parser.add_argument(
        "-o", "--output_gtf_file",
        dest="output_gtf_f",
        required=True,
        help="Output GTF file"
    )

    parser.add_argument(
        "-r", "--ref_gtf_file",
        dest="ref_gtf_f",
        required=False,
        default=None,
        help="Reference GTF file"
    )
    parser.add_argument(
        "-R", "--ref_db_file",
        dest="ref_db_f",
        required=False,
        default=None,
        help="Reference DB file"
    )

    parser.add_argument(
        "-S", "--select_transcripts_f",
        dest="select_transcripts_f",
        required=False,
        default=None,
        help="File containing selected reference transcript IDs"
    )

    parser.add_argument(
        "-s", "--select_transcripts",
        dest="select_transcripts_list",
        required=False,
        default=None,
        nargs='+',
        help="List of selected reference transcript IDs"
    )



    args = parser.parse_args()

    print("Output GTF file:", args.output_gtf_f)

    if args.query_gtf_fs is not None:
        print("Query GTF files:",','.join(args.query_gtf_fs) )
    if args.query_gtf_fofn is not None:
        print("Query FOFN:",args.query_gtf_fofn)
    if args.ref_gtf_f is not None:
        print("Reference GTF file:", args.ref_gtf_f)
    if args.ref_db_f is not None:
        print("Reference DB file:", args.ref_db_f)
    if args.select_transcripts_f is not None:
        print("Select reference transcripts file:", args.select_transcripts_f)
    if args.select_transcripts_list is not None:
        print("Select reference transcripts list:", ', '.join(args.select_transcripts_list))


    query_gtf_fs=args.query_gtf_fs
    query_gtf_fofn=args.query_gtf_fofn
    output_gtf_f=args.output_gtf_f
    ref_gtf_f=args.ref_gtf_f
    ref_db_f=args.ref_db_f
    select_transcripts_f=args.select_transcripts_f
    select_transcripts_list=args.select_transcripts_list




    #Composite arguments
    ref_given=(ref_gtf_f is not None) or (ref_db_f is not None)
    subset_given=(select_transcripts_f is not None) or (select_transcripts_list is not None)
    query_given=(query_gtf_fs is not None) or (query_gtf_fofn is not None)
    ######
    #Based on the combination of arguments, you should use different steps. Here as a prototype, we do a full union of all databases (ref+query)
    #####
    if not query_given:
        raise ValueError('Query not given either as list or FOFN')
    if query_given and (query_gtf_fs is None):
        query_gtf_fs=pd.read_csv(query_gtf_fofn,sep='\t',header=None)
        query_gtf_fs=query_gtf_fs[0].tolist()


    process_ref_db(ref_db_f,ref_gtf_f)
    query_dbs=create_feature_db_from_gtf(query_gtf_fs)
    #Based on the combination of arguments, you should use different steps. Here as a prototype, we do a full union of all databases (ref+query)
    if ref_given and not subset_given:
        query_dbs=process_ref_db(ref_db_f,ref_gtf_f)+query_dbs

    #Merging, formatting and writing. Shouldn't depend on strategy and only base on gene_dict
    merged_gene_dict=merge_dbs(query_dbs)
    output_str=format_gtf_output(merged_gene_dict)
    with open(output_gtf_f, "w") as out_f:
        out_f.write(output_str)


if __name__=='__main__':
    main()


#USAGE: python collect_gtfs.py -q /lustre/scratch126/humgen/projects/isogut/workdirs/isoquant_twopass/8e/0438868efe618be3643e238ebba8ec/chr12_99988580_106538194/chr12_99988580_106538194_renamed/chr12_99988580_106538194.transcript_models.gtf /lustre/scratch126/humgen/projects/isogut/workdirs/isoquant_twopass/42/df19e313a1f8128fab2349f5d28005/chr8_129706194_142880609/chr8_129706194_142880609_renamed/chr8_129706194_142880609.transcript_models.gtf -R /lustre/scratch126/humgen/projects/isogut/utils/gencode.v46.annotation.sorted.gtf.db  -o merge.gtf
