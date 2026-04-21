include {GT_MATCH_POOL_AGAINST_PANEL; ASSIGN_DONOR_FROM_PANEL; ASSIGN_DONOR_OVERALL; VIREO_GT_FIX_HEADER; COMBINE_ASSIGN} from  '../../modules/match_gt.nf'
include {MATCH_BARCODES; COMBINE_DONOR_ASIGNMENTS; REMOVE_DONOR_ASIGNMENTS; REPLACE_DONOR_ASIGNMENTS} from '../../modules/matching_barcodes.nf'

workflow DONOR_ASSIGNMENT {
  take:
    vireo_map
    sample_donor_vcf_ch
  main:
    Channel.empty().set { ch_versions }
//gt and barcode matching
    vireo_map2=vireo_map
    if (params.run_gtcheck == 'TRUE'){
      VIREO_GT_FIX_HEADER(sample_donor_vcf_ch, params.genome_fasta_f)
      sample_gt_vcf_ch=VIREO_GT_FIX_HEADER.out.gt_pool

      ref_vcf_ch = Channel.fromPath(
        params.tsv_donor_panel_vcfs,
        followLinks: true,
        checkIfExists: true
      ).splitCsv(header: true, sep: '\t')
      .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }

      def tsv_entries = file(params.tsv_donor_panel_vcfs).splitCsv(header: true, sep: '\t')
      def has_multiple_ref_vcfs = tsv_entries.size() > 1

      gt_math_pool_against_panel_input=sample_gt_vcf_ch.combine(ref_vcf_ch)
        .map { pool_id, vcf, vcf_tbi, ref_label, ref_vcf, ref_csi -> 
          tuple(pool_id, vcf, vcf_tbi, ref_label, ref_vcf, ref_csi)
        }

      // now match genotypes against a panels
      GT_MATCH_POOL_AGAINST_PANEL(gt_math_pool_against_panel_input)
      ch_versions = ch_versions.mix(GT_MATCH_POOL_AGAINST_PANEL.out.versions)

      // group by panel id
      GT_MATCH_POOL_AGAINST_PANEL.out.gtcheck_results.unique()
        .groupTuple()
        .set { gt_check_by_panel }

      ASSIGN_DONOR_FROM_PANEL(gt_check_by_panel)
      ch_versions = ch_versions.mix(ASSIGN_DONOR_FROM_PANEL.out.versions)
      ASSIGN_DONOR_FROM_PANEL.out.gtcheck_assignments.unique()
        .groupTuple()
        .set{ ch_donor_assign_panel }

      assignment_ch=ASSIGN_DONOR_FROM_PANEL.out.gtcheck_assignments

      ASSIGN_DONOR_OVERALL(ch_donor_assign_panel)
      ch_versions = ch_versions.mix(ASSIGN_DONOR_OVERALL.out.versions)
      assignment_ch=assignment_ch.mix(ASSIGN_DONOR_OVERALL.out.donor_match_table_with_pool_id)

      COMBINE_ASSIGN(assignment_ch.map { val, path -> path }.collect())
      ch_versions = ch_versions.mix(COMBINE_ASSIGN.out.versions)
      match_gt=COMBINE_ASSIGN.out.gt_match_table

      vireo_map2=vireo_map2.mix(COMBINE_ASSIGN.out.gt_match_table)
    }
    if (params.run_barcode_check == 'TRUE'){ 
      vireo_map=vireo_map.map { file -> tuple('vireo_map', file) }
      MATCH_BARCODES(vireo_map, params.ref_barcode_list, params.min_n_common_barcodes, params.min_common_barcodes_percent)
      vireo_map2=vireo_map2.mix(MATCH_BARCODES.out.match_barcodes)
    }
    vireo_map2 = vireo_map2.collect()

    COMBINE_DONOR_ASIGNMENTS(vireo_map2)
}