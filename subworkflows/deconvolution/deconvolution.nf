include {DONOR_ASSIGNMENT} from '../donor_assignment/donor_assignment.nf'
include { MPILEUP; CELLSNP; VIREO; SUBSET_VCF }  from '../../modules/deconvolution.nf'
//include {GET_BARCODES; SPLIT_BAM; INDEX_SPLIT_BAM} from '../../modules/split_bam.nf'
include {SPLIT_BAM_SINTO; INDEX_SPLIT_BAM} from '../../modules/split_bam.nf'

workflow DECONVOLUTION {
  take:
    mapped_reads
  main:
    Channel.empty().set { ch_versions }
    if (mapped_reads == 'independent workflow') {
    Channel
      .fromPath(params.input_samples_path)
      .splitCsv(sep: ',', header: true)
      //.filter { it -> !(it.sample_id in excluded_samples_list) }
      .map { it ->
        def sample_id = it.sample_id
        def bam_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam"
        def bai_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam.bai"
        [sample_id, bam_path, bai_path]
      }
      .set { mapped_reads }
    }

    Channel
      .fromPath(params.input_samples_path)
      .splitCsv(sep: ',', header: true)
      .map { it -> [it.sample_id, it.nr_samples_multiplexed] } // Create a tuple with bam_path and sample_id
      .set { sampleNames_nrDons }

    // Only perform deconvolution on samples that are more than 1 donor.
    combined_ch_deconv =  mapped_reads.combine(sampleNames_nrDons, by: 0)
    combined_ch_deconv.filter { experiment, bam, bai, npooled -> npooled == '1' }.map { experiment, bam, bai, _ -> [experiment, bam, bai] }.set{not_for_deconv}
    combined_ch_deconv.filter { experiment, bam, bai, npooled -> npooled != '1' }.map { experiment, bam, bai, _ -> [experiment, bam, bai] }.set{for_deconv}

    MPILEUP(for_deconv, params.genome_fasta_f)
    SUBSET_VCF(MPILEUP.out, params.subset_regions_bed)
    CELLSNP(SUBSET_VCF.out)

    sampleNames_cellsnp_nrDons = CELLSNP.out.combine(sampleNames_nrDons, by: 0)
    VIREO(sampleNames_cellsnp_nrDons)
    mapped_reads = mapped_reads.combine(VIREO.out.barcode_list, by: 0)
    SPLIT_BAM_SINTO(mapped_reads)
    donor_bams = SPLIT_BAM_SINTO.out.donor_bams.flatMap { val, paths ->
      paths.collect { p -> tuple(p.simpleName, p) } }
    INDEX_SPLIT_BAM(donor_bams)
    fullBam_ch_pre=INDEX_SPLIT_BAM.out.donor_bams
    fullBam_ch=fullBam_ch_pre.mix(not_for_deconv)
    if (params.run_gtcheck== 'TRUE' || params.run_barcode_check == 'TRUE') {
        VIREO.out.sample_donor_ids.collectFile(name: 'vireo_map.tsv', newLine: true) { sample_id, file_path ->
            "${sample_id}\t${file_path}"
        }
        .set { vireo_map }
        DONOR_ASSIGNMENT(vireo_map, VIREO.out.sample_donor_vcf)
    }

  emit:
    fullBam_ch
}
