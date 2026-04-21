include { SPLIT_READS; REMOVE_PRIMER; TAG_BAM; REFINE_READS } from '../../modules/fltnc.nf'
include {BARCODE_CORRECTION; GET_BARCODES; SUPSET_BAM; DEDUP_READS; COMBINE_DEDUPS; COMBINE_MUPPED; BAM_STATS} from '../../modules/barcodes.nf'

include { PBMM2 } from '../../modules/pbmm2.nf'

workflow BAM_PROCESSING {
    main:
      /// Obtaining (sample_id,bam_file) tuples from the input_samples.csv file
      Channel
        .fromPath(params.input_samples_path)
        .splitCsv(sep: ',', header: true)
        .map { it -> [it.sample_id, it.long_read_path] } // Create a tuple with bam_path and sample_id
        .set { hifi_bam_tuples }

      /// Every process from now on outputs a (sample_id,bam_file) tuple which is fed on to the next process
      SPLIT_READS(hifi_bam_tuples, params.skera_primers)
      REMOVE_PRIMER(SPLIT_READS.out.split_reads_tuple, params.tenx_primers)
      TAG_BAM(REMOVE_PRIMER.out.removed_primer_tuple)
      REFINE_READS(TAG_BAM.out.tagged_tuple, params.tenx_primers, params.min_polya_length)
      //barcode correction
      BARCODE_CORRECTION(REFINE_READS.out.refined_reads, params.threeprime_whitelist, params.barcode_correction_method, params.barcode_correction_percentile)
      GET_BARCODES(BARCODE_CORRECTION.out.barcode_corrected_tuple, params.number_of_chunks)
      barcode_channel=GET_BARCODES.out.barcodes_tuple.transpose()
      combined_ch = BARCODE_CORRECTION.out.barcode_corrected_tuple.combine(barcode_channel, by: 0)
      SUPSET_BAM(combined_ch)
      DEDUP_READS(SUPSET_BAM.out.chunk_tuple, params.dedup_batch_size)
      //these three should create combined dedup files and their stats
      deduped_chunks=DEDUP_READS.out.dedup_tuple
      deduped_chunks_ch=deduped_chunks.groupTuple()
      COMBINE_DEDUPS(deduped_chunks_ch)
      BAM_STATS(COMBINE_DEDUPS.out.dedup_tuple, params.barcode_correction_method, params.barcode_correction_percentile)
      //mapping
      PBMM2(DEDUP_READS.out.dedup_tuple, params.genome_fasta_f)
      mapped_chunks_ch=PBMM2.out.map_tuple.groupTuple()
      COMBINE_MUPPED(mapped_chunks_ch)
    emit:
    mapped_reads = COMBINE_MUPPED.out.combined_bam_tuple
}
