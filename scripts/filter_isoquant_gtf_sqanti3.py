import pandas as pd
import argparse

# Initialize argument parser
parser = argparse.ArgumentParser(description="Filter a GTF file by associated genes and transcripts.")
parser.add_argument("--gtf", required=True, help="Path to the input GTF file.")
parser.add_argument("--filter", required=True, help="Path to the filter list CSV file.")
parser.add_argument("--output", required=True, help="Path to save the filtered GTF file.")
args = parser.parse_args()

# Load the filter list
filter_list = pd.read_csv(args.filter,sep='\t')
filter_list=filter_list.loc[filter_list['filter_result']=='Isoform',['associated_gene','associated_transcript']]

# Load the GTF file
gtf_columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gtf = pd.read_csv(args.gtf, sep="\t", names=gtf_columns, comment="#", header=None)

# Function to extract gene_id and transcript_id from the attribute field
def extract_ids(attribute):
    gene_id = transcript_id = None
    for entry in attribute.split(";"):
        if "gene_id" in entry:
            gene_id = entry.split('"')[1]
        elif "transcript_id" in entry:
            transcript_id = entry.split('"')[1]
    return gene_id, transcript_id

# Apply extraction to the GTF DataFrame
gtf[["gene_id", "transcript_id"]] = gtf["attribute"].apply(lambda x: pd.Series(extract_ids(x)))

# Separate gene entries
gene_entries = gtf[gtf["feature"] == "gene"]

# Filter gene entries by associated_gene
filtered_genes = gene_entries[gene_entries["gene_id"].isin(filter_list["associated_gene"])]

# Filter non-gene entries by associated_gene and associated_transcript
non_gene_entries = gtf[gtf["feature"] != "gene"]
filtered_non_genes = non_gene_entries.merge(filter_list, left_on=["gene_id", "transcript_id"], right_on=["associated_gene", "associated_transcript"])

# Count remaining transcripts per gene
remaining_transcripts = filtered_non_genes["gene_id"].value_counts()

# Update the "transcripts" attribute in gene records #TODO: this doesn't work to update the number of gene transcripts. FIX
def update_transcripts(attribute, gene_id):
    if "transcripts" in attribute:
        count = remaining_transcripts.get(gene_id, 0)
        print(gene_id)
        print(count)
        print('-------')
        attribute = attribute.replace(f'transcripts "{attribute.split("transcripts ")[1].split(";")[0]}"', f'transcripts "{count}"')
    return attribute

#filtered_genes["attribute"] = filtered_genes.apply(lambda row: update_transcripts(row["attribute"], row["gene_id"]), axis=1)

# Combine the filtered DataFrames and preserve original order
filtered_gtf = pd.concat([filtered_genes, filtered_non_genes]).sort_index()

# Drop temporary columns and save the filtered GTF file
filtered_gtf.drop(columns=["gene_id", "transcript_id","associated_gene","associated_transcript"], inplace=True)
filtered_gtf.to_csv(args.output, sep="\t", index=False, header=False, quoting=3, quotechar='', escapechar='')
print(filtered_gtf)
print(f"Filtered GTF saved as '{args.output}'")
