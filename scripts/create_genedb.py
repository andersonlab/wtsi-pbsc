import gffutils
import argparse

parser = argparse.ArgumentParser(description='Converts GFF to gffutils database used by isoquant')

parser.add_argument('-g', '--gff',dest='gff_f',required=True,type=str)
parser.add_argument('-o', '--db',dest='output_f',required=True,type=str)
args = parser.parse_args()

gff_f=args.gff_f
output_f=args.output_f
db = gffutils.create_db(gff_f, dbfn=output_f, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True,disable_infer_genes=True,disable_infer_transcripts=True)
