

#USAGE: python /nfs/team152/oe2/isogut/scripts/workflows/isoseq/scripts/gene_coords_by_biotype.py -d /lustre/scratch126/humgen/projects/isogut/utils/ensembl/Homo_sapiens.GRCh38.113.gtf.db -b IG_C_gene IG_D_gene IG_J_gene IG_LV_gene IG_V_gene TR_C_gene TR_J_gene TR_V_gene TR_D_gene -o vdj_genes.bed
import argparse
import pybedtools
import gffutils
import os

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Extract and merge genes from a GTF or GFFutils DB file based on biotypes.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--gtf", required=False, help="Path to the Ensembl GTF file (if no DB file is provided).")
    group.add_argument("-d", "--db" , required=False, help="Path to an existing GFFutils database file (optional alternative to GTF).")

    parser.add_argument("-b", "--biotypes",nargs='+', required=True, help="Gene biotypes to filter.")
    parser.add_argument("-o", "--output", required=True, help="Output BED file for merged gene regions.")

    return parser.parse_args()

def load_or_create_db(gtf_file, db_file):
    """Load an existing GFFutils database or create one from a GTF file."""
    if db_file and os.path.exists(db_file):
        print(f"Loading existing database: {db_file}")
        return gffutils.FeatureDB(db_file, keep_order=True)
    elif gtf_file:
        print(f"Creating database from GTF: {gtf_file}")
        return gffutils.create_db(gtf_file, dbfn=":memory:", force=True, keep_order=True,
                                  disable_infer_transcripts=True, disable_infer_genes=True)
    else:
        raise ValueError("Either a GTF file or a GFFutils DB file must be provided.")

def main():
    args = parse_args()
    biotypes = set(args.biotypes)

    # Load existing database or create from GTF
    db = load_or_create_db(args.gtf, args.db)

    # Extract gene coordinates for selected biotypes
    gene_coords = []
    for gene in db.features_of_type("gene"):
        if "gene_biotype" in gene.attributes:
            biotype = gene.attributes["gene_biotype"][0]  # Extract biotype
            if biotype in biotypes:
                chrom=gene.chrom if gene.chrom.startswith('chr') else 'chr'+gene.chrom
                gene_coords.append([chrom, str(gene.start - 1), str(gene.end)])  # Convert to BED format (0-based start)

    if not gene_coords:
        print("No matching genes found for the specified biotypes.")
        return

    # Convert to BED format and process with pybedtools
    bed_lines = ['\t'.join(coord) for coord in gene_coords]
    bedtool = pybedtools.BedTool('\n'.join(bed_lines), from_string=True)

    # Merge overlapping/adjacent regions
    merged = bedtool.sort().merge()

    # Save to output file
    merged.saveas(args.output)
    print(f"Saved merged BED file to {args.output}")

if __name__ == "__main__":
    main()
