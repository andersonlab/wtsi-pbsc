include {GFFREAD_PROTEINS; HMMSCAN_PFAM} from '../../modules/pfam_annotation.nf'

workflow PFAM_ANNOTATION_WF {
  take:
    filtered_gtf

  main:
    if (filtered_gtf == 'independent workflow') {
      Channel
        .fromPath("${params.results_output}results/transcript_info/sqanti3/sqanti3_filter/*.filtered.gtf")
        .set { filtered_gtf }
    }
    GFFREAD_PROTEINS(filtered_gtf, params.genome_fasta_f)
    HMMSCAN_PFAM(GFFREAD_PROTEINS.out.proteins_faa)
}
