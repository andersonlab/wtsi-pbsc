import argparse
import logging
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

PROGRESS_EVERY = 1000

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
TX_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def count_gtf_features(path):
    counts = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            feature = line.split('\t', 3)[2]
            counts[feature] = counts.get(feature, 0) + 1
    return counts


def process_file(path, seen_gene_ids, seen_tx_ids, genes, counters):
    current_tx = None
    tx_is_dup = True

    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            feature = fields[2]

            if feature == 'gene':
                gene_id = GENE_ID_RE.search(line).group(1)
                if gene_id in seen_gene_ids:
                    counters['dup_genes'] += 1
                    continue
                seen_gene_ids.add(gene_id)
                genes[gene_id] = {
                    'gene_line': line,
                    'chrom': fields[0],
                    'start': int(fields[3]),
                    'transcripts': [],
                }
                counters['genes'] += 1
                continue

            if feature == 'transcript':
                tx_id = TX_ID_RE.search(line).group(1)
                tx_is_dup = tx_id in seen_tx_ids
                if tx_is_dup:
                    counters['dup_transcripts'] += 1
                    continue
                seen_tx_ids.add(tx_id)

                gene_id = GENE_ID_RE.search(line).group(1)
                gene = genes.get(gene_id)
                if gene is None:
                    # Defensive fallback -- shouldn't happen with well-formed GENCODE/IsoQuant
                    # input, where every transcript's gene line always precedes it.
                    logging.warning(
                        'transcript_id "%s" has no preceding gene_id "%s" record; synthesizing one',
                        tx_id, gene_id,
                    )
                    gene_line = '\t'.join([
                        fields[0], fields[1], 'gene', fields[3], fields[4], '.', fields[6], '.',
                        'gene_id "{0}";\n'.format(gene_id),
                    ])
                    gene = {'gene_line': gene_line, 'chrom': fields[0], 'start': int(fields[3]), 'transcripts': []}
                    genes[gene_id] = gene
                    seen_gene_ids.add(gene_id)
                    counters['genes'] += 1

                current_tx = {'start': int(fields[3]), 'lines': [line]}
                gene['transcripts'].append(current_tx)
                counters['transcripts'] += 1
                if counters['transcripts'] % PROGRESS_EVERY == 0:
                    logging.info('Kept %d transcripts so far..', counters['transcripts'])
                continue

            # Any other feature row (exon, CDS, UTR, start_codon, stop_codon, ...) belongs to
            # whichever transcript line preceded it -- copied through verbatim, untouched.
            if not tx_is_dup:
                current_tx['lines'].append(line)


def main():
    parser = argparse.ArgumentParser(
        description='Concatenate a reference GTF with one or more query GTFs, keeping the '
                    'first-seen record for any gene_id/transcript_id that appears in more than '
                    'one source (reference wins, then query files in the given order), and '
                    'write the result in GENCODE-style order: genes sorted by (chrom, start), '
                    'each gene\'s transcripts (sorted by start) kept fully contiguous with it. '
                    'Every kept line is copied verbatim from its source -- no reformatting, no '
                    'dropped feature types.'
    )
    parser.add_argument('-r', '--ref-gtf', dest='ref_gtf_f', required=True, help='Reference GTF (highest priority: wins any gene_id/transcript_id conflict)')
    parser.add_argument('-q', '--query-gtfs', dest='query_gtf_fs', required=True, nargs='+', help='Query GTFs, in priority order after the reference')
    parser.add_argument('-o', '--out', dest='out_f', required=True, help='Output GTF path')
    args = parser.parse_args()

    seen_gene_ids = set()
    seen_tx_ids = set()
    genes = {}
    counters = {'genes': 0, 'transcripts': 0, 'dup_genes': 0, 'dup_transcripts': 0}

    process_file(args.ref_gtf_f, seen_gene_ids, seen_tx_ids, genes, counters)
    for query_gtf_f in args.query_gtf_fs:
        process_file(query_gtf_f, seen_gene_ids, seen_tx_ids, genes, counters)

    logging.info(
        'Kept %d genes (%d duplicate gene_id records dropped) and %d transcripts '
        '(%d duplicate transcript_id records dropped)',
        counters['genes'], counters['dup_genes'], counters['transcripts'], counters['dup_transcripts'],
    )

    with open(args.out_f, 'w') as out_f:
        for gene_id, gene in sorted(genes.items(), key=lambda kv: (kv[1]['chrom'], kv[1]['start'])):
            out_f.write(gene['gene_line'])
            for tx in sorted(gene['transcripts'], key=lambda t: t['start']):
                out_f.writelines(tx['lines'])

    logging.info('Wrote %s', args.out_f)

    # Verification: independently re-read the written file and confirm its feature counts match
    # what we intended to write (len(seen_*_ids)) -- catches write-time issues (e.g. a bucket
    # dropped during the sort/write pass), not just in-memory bookkeeping.
    feature_counts = count_gtf_features(args.out_f)
    actual_gene_count = feature_counts.get('gene', 0)
    actual_tx_count = feature_counts.get('transcript', 0)
    expected_gene_count = len(seen_gene_ids)
    expected_tx_count = len(seen_tx_ids)

    logging.info(
        f"Verification: expected {expected_gene_count} genes, found {actual_gene_count} in {args.out_f}"
    )
    logging.info(
        f"Verification: expected {expected_tx_count} transcripts, found {actual_tx_count} in {args.out_f}"
    )

    if actual_gene_count != expected_gene_count or actual_tx_count != expected_tx_count:
        raise RuntimeError(
            f"{args.out_f} feature counts do not match expectations: "
            f"genes expected={expected_gene_count} actual={actual_gene_count}, "
            f"transcripts expected={expected_tx_count} actual={actual_tx_count}"
        )
    logging.info("Verification PASSED: extended_annotation.gtf feature counts match expectations")


if __name__ == '__main__':
    main()
