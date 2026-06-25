include {run_isoquant_firstPass; run_isoquant_firstPass_perSample; create_model_construction_bam_perChr; run_isoquant_chunked_merged; replace_novel_names; collect_gtfs} from '../../modules/isoquant.nf'
include {chroms} from '../core/chroms.nf'
include {run_isoquant_firstPass_withmodelconstruction; replace_novel_names_firsPass_singlenovelname} from '../../modules/isoquant.nf'
include {collect_counts_as_mtx_perChr as collect_isoform_counts_as_mtx_perChr} from '../../modules/isoquant.nf'
include {collect_counts_as_mtx_perChr as collect_gene_counts_as_mtx_perChr} from '../../modules/isoquant.nf'
include {collect_counts_as_mtx_perChr as collect_intron_include_counts_as_mtx_perChr} from '../../modules/isoquant.nf'
include {collect_counts_as_mtx_perChr as collect_intron_exclude_counts_as_mtx_perChr} from '../../modules/isoquant.nf'
include {collect_counts_as_mtx_perChr as collect_exon_include_counts_as_mtx_perChr} from '../../modules/isoquant.nf'
include {collect_counts_as_mtx_perChr as collect_exon_exclude_counts_as_mtx_perChr} from '../../modules/isoquant.nf'

include {collect_counts_as_mtx_perChr as collect_isoform_chrM_counts_as_mtx_perChr} from '../../modules/isoquant.nf'
include {collect_counts_as_mtx_perChr as collect_gene_chrM_counts_as_mtx_perChr} from '../../modules/isoquant.nf'

include {collect_mtx_as_h5ad as collect_isoform_mtx_as_h5ad} from '../../modules/isoquant.nf'
include {collect_mtx_as_h5ad as collect_gene_mtx_as_h5ad} from '../../modules/isoquant.nf'
include {collect_mtx_as_h5ad as collect_intron_include_mtx_as_h5ad} from '../../modules/isoquant.nf'
include {collect_mtx_as_h5ad as collect_intron_exclude_mtx_as_h5ad} from '../../modules/isoquant.nf'
include {collect_mtx_as_h5ad as collect_exon_include_mtx_as_h5ad} from '../../modules/isoquant.nf'
include {collect_mtx_as_h5ad as collect_exon_exclude_mtx_as_h5ad} from '../../modules/isoquant.nf'

include {collect_mtx_as_h5ad as collect_isoform_chrM_mtx_as_h5ad} from '../../modules/isoquant.nf'
include {collect_mtx_as_h5ad as collect_gene_chrM_mtx_as_h5ad} from '../../modules/isoquant.nf'

include {format_intron_exon_grouped_counts_perChr as format_intron_grouped_counts_firstPass_perChr} from '../../modules/isoquant.nf'
include {format_intron_exon_grouped_counts_perChr as format_intron_grouped_counts_secondPass_perChr} from '../../modules/isoquant.nf'
include {format_intron_exon_grouped_counts_perChr as format_exon_grouped_counts_firstPass_perChr} from '../../modules/isoquant.nf'
include {format_intron_exon_grouped_counts_perChr as format_exon_grouped_counts_secondPass_perChr} from '../../modules/isoquant.nf'

include {find_mapped_and_unmapped_regions_perChr; suggest_splits_binarySearch_merged; split_bam_perChunk} from '../../modules/smartSplit.nf'

workflow isoquant_twopass_perChr_wf {
  take:
    isoquant_preprocess_bam_perChr_ch
  main:
    isoquant_firstpass_output_ch=run_isoquant_firstPass(isoquant_preprocess_bam_perChr_ch, params.genedb_f, params.genome_fasta_f, params.genome_fasta_f + ".fai")
    isoquant_firstpass_output_ch
    .map{ chrom,sample_id,isoquant_output_dir,read_assignment_f,bam -> [chrom,sample_id,read_assignment_f,bam] }
    .set{ model_construction_bam_input_ch }
    model_construction_bam_ch=create_model_construction_bam(model_construction_bam_input_ch)

    model_construction_bam_ch
    .groupTuple(by:0)
    .set{ isoquant_secondpass_input_ch }

    isoquant_secondpass_output_ch=run_isoquant_perChr(isoquant_secondpass_input_ch, params.genedb_f, params.genome_fasta_f)

  emit:
    isoquant_secondpass_output_ch
}



workflow collect_exon_intron_coutns_perChr_wf {
  take:
    isoquant_firstpass_output_ch
    isoquant_output_novel_names_ch
    bam_nums_perChr_ch
    chunks
  main:
    /////////////////////////////////////////////////////////
    //////COLLECTION INTRON/EXON GROUPED COUNTS////////////
    /////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////
    //////A)Collecting intron grouped counts as MTX//////////
    /////////////////////////////////////////////////////////


    isoquant_firstpass_output_ch
    ///.filter{tpl -> (tpl[0]=='chr2') && (tpl[1]=='Isogut15045390') }
    .map{chrom,sample_id,isoquant_output_dir,read_assignments_f,processed_bam -> [chrom,"${sample_id}.${chrom}.intron_grouped_counts","${isoquant_output_dir}/${sample_id}.${chrom}/${sample_id}.${chrom}.intron_grouped_tag_CB_counts.tsv"]} | format_intron_grouped_counts_firstPass_perChr | set{firstPass_intron_grouped_formatted_ch}

    isoquant_output_novel_names_ch
    ///.filter{tpl -> tpl[1]=='chr2_137822130_165323855'}
    .map{chrom,programmaticRegion,isoquant_output_dir -> [chrom,"${programmaticRegion}.intron_grouped_counts","${isoquant_output_dir}/${programmaticRegion}.intron_grouped_tag_CB_counts.tsv"]} | format_intron_grouped_counts_secondPass_perChr | set{secondPass_intron_grouped_formatted_ch}

    firstPass_intron_grouped_formatted_ch
        .combine(bam_nums_perChr_ch,by:0)
        .map{chrom, include_f,exclude_f,chrom_sample_size -> [groupKey(chrom,chrom_sample_size),include_f]}
        .groupTuple(by:0)
        .map {tpl -> [tpl[0],tpl[1]]}
        .combine(
        secondPass_intron_grouped_formatted_ch
            .map{chrom, prefix,intron_f -> [groupKey(chrom,chunks),intron_f]}
            .groupTuple(by:0)
            .map {tpl -> [tpl[0],tpl[1]]}
        ,by:0
        )
        .map { chrom, firstPasslist, secondPasslist -> [chrom, firstPasslist + secondPasslist] } | collect_intron_include_counts_as_mtx_perChr | set{intron_include_mtx}

        collect_intron_include_mtx_as_h5ad(intron_include_mtx.chrom_mtx | collect, 'introns_include')

    firstPass_intron_grouped_formatted_ch
        .combine(bam_nums_perChr_ch,by:0)
        .map{chrom, include_f,exclude_f,chrom_sample_size -> [groupKey(chrom,chrom_sample_size),exclude_f]}
        .groupTuple(by:0)
        .map {tpl -> [tpl[0],tpl[1]]}
        .combine(
        secondPass_intron_grouped_formatted_ch
            .map{chrom, prefix,intron_f -> [groupKey(chrom,chunks),intron_f]}
            .groupTuple(by:0)
            .map {tpl -> [tpl[0],tpl[1]]}
        ,by:0
        )
        .map { chrom, firstPasslist, secondPasslist -> [chrom, firstPasslist + secondPasslist] } | collect_intron_exclude_counts_as_mtx_perChr | set{intron_exclude_mtx}

        collect_intron_exclude_mtx_as_h5ad(intron_exclude_mtx.chrom_mtx | collect, 'introns_exclude')




        /////////////////////////////////////////////////////////
        //////A)Collecting exon grouped counts as MTX//////////
        /////////////////////////////////////////////////////////

        isoquant_firstpass_output_ch
        ///.filter{tpl -> (tpl[0]=='chr2') && (tpl[1]=='Isogut15045390') }
        .map{chrom,sample_id,isoquant_output_dir,read_assignments_f,processed_bam -> [chrom,"${sample_id}.${chrom}.exon_grouped_counts","${isoquant_output_dir}/${sample_id}.${chrom}/${sample_id}.${chrom}.exon_grouped_tag_CB_counts.tsv"]} | format_exon_grouped_counts_firstPass_perChr | set{firstPass_exon_grouped_formatted_ch}

        isoquant_output_novel_names_ch
        ///.filter{tpl -> tpl[1]=='chr2_137822130_165323855'}
        .map{chrom,programmaticRegion,isoquant_output_dir -> [chrom,"${programmaticRegion}.exon_grouped_counts","${isoquant_output_dir}/${programmaticRegion}.exon_grouped_tag_CB_counts.tsv"]} | format_exon_grouped_counts_secondPass_perChr | set{secondPass_exon_grouped_formatted_ch}


        firstPass_exon_grouped_formatted_ch
            .combine(bam_nums_perChr_ch,by:0)
            .map{chrom, include_f,exclude_f,chrom_sample_size -> [groupKey(chrom,chrom_sample_size),include_f]}
            .groupTuple(by:0)
            .map {tpl -> [tpl[0],tpl[1]]}
            .combine(
            secondPass_exon_grouped_formatted_ch
                .map{chrom, prefix,exon_f -> [groupKey(chrom,chunks),exon_f]}
                .groupTuple(by:0)
                .map {tpl -> [tpl[0],tpl[1]]}
            ,by:0
            )
            .map { chrom, firstPasslist, secondPasslist -> [chrom, firstPasslist + secondPasslist] } | collect_exon_include_counts_as_mtx_perChr | set{exon_include_mtx}

            collect_exon_include_mtx_as_h5ad(exon_include_mtx.chrom_mtx | collect, 'exons_include_mtx')

        firstPass_exon_grouped_formatted_ch
            .combine(bam_nums_perChr_ch,by:0)
            .map{chrom, include_f,exclude_f,chrom_sample_size -> [groupKey(chrom,chrom_sample_size),exclude_f]}
            .groupTuple(by:0)
            .map {tpl -> [tpl[0],tpl[1]]}
            .combine(
            secondPass_exon_grouped_formatted_ch
                .map{chrom, prefix,exon_f -> [groupKey(chrom,chunks),exon_f]}
                .groupTuple(by:0)
                .map {tpl -> [tpl[0],tpl[1]]}
            ,by:0
            )
            .map { chrom, firstPasslist, secondPasslist -> [chrom, firstPasslist + secondPasslist] } | collect_exon_exclude_counts_as_mtx_perChr | set{exon_exclude_mtx}

            collect_exon_exclude_mtx_as_h5ad(exon_exclude_mtx.chrom_mtx | collect, 'exons_exclude_mtx')

    /////////////////////////////////////////////////////////
    //////END: COLLECTION INTRON/EXON GROUPED COUNTS////////////
    /////////////////////////////////////////////////////////
}


////////////////////////////////
//////chrM workflows///////////
///////////////////////////////
workflow isoquant_chrM {
  take:
    isoquant_preprocess_bam_perChr_ch
    chrom_sizes_f
  main:
    ///////////////////////////////////////////////////////
    //////////////////A-FIRST PASS/////////////////////////
    ///////////////////////////////////////////////////////
    isoquant_preprocess_bam_perChr_ch.groupTuple(by:0).map{tpl -> [tpl[0],(tpl[1]).size()]}.set{bam_nums_perChr_ch}
    isoquant_firstpass_output_ch=run_isoquant_firstPass_withmodelconstruction(isoquant_preprocess_bam_perChr_ch, params.genedb_f, params.genome_fasta_f, params.genome_fasta_f + ".fai")

    isoquant_firstpass_output_ch.map{chrom,sample_id,isoquant_tar,isoquant_input_bam -> [chrom,sample_id,isoquant_tar]}
    .set{replace_novel_names_input_ch}


    //Updating names of novel transcript so they don't clash between chunks
    isoquant_output_novel_names_ch=replace_novel_names_firsPass_singlenovelname(replace_novel_names_input_ch)

    ///////////////////////////////////////////////////////
    //////////////////END: FIRST PASS//////////////////////
    ///////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    //////////////////B-OUTPUT CHANNELs/////////////////////////
    ////////////////////////////////////////////////////////////
    // Each count channel emits [chrom, [tars], [patterns], [needs_chr_filter]].
    // chrM tars are always per-chromosome (no filtering needed), so needs_chr_filter is always false.
    isoquant_output_novel_names_ch
      .combine(bam_nums_perChr_ch,by:0)
      .map{chrom,sample_id,tar,n -> [groupKey(chrom,n),tar,"*.discovered_transcript_grouped_tag_CB_counts.linear.tsv",false]}
      .groupTuple(by:0)
      .map{chrom,tars,patterns,filters -> [chrom,tars,patterns,filters]}
      .set{output_isoform_counts_ch}

    isoquant_output_novel_names_ch
      .combine(bam_nums_perChr_ch,by:0)
      .map{chrom,sample_id,tar,n -> [groupKey(chrom,n),tar,"*.discovered_gene_grouped_tag_CB_counts.linear.tsv",false]}
      .groupTuple(by:0)
      .map{chrom,tars,patterns,filters -> [chrom,tars,patterns,filters]}
      .set{output_gene_counts_ch}

    isoquant_output_novel_names_ch
      .combine(bam_nums_perChr_ch,by:0)
      .map{chrom,sample_id,tar,n -> [groupKey(chrom,n),tar]}
      .groupTuple(by:0)
      .set{output_existing_gtf_ch}

    output_corrected_reads_ch = Channel.empty()
    output_assignment_reads_ch = Channel.empty()
    output_transcriptmodel_reads_ch = Channel.empty()


    ///////////////////////////////////////////////////////
    //////////////////END: OUTPUT CHANNELS/////////////////
    ///////////////////////////////////////////////////////

  emit:
    isoform_counts=output_isoform_counts_ch
    gene_counts=output_gene_counts_ch
    existing_gtf=output_existing_gtf_ch
    assignment_reads=output_assignment_reads_ch
    transcriptmodel_reads=output_transcriptmodel_reads_ch
    corrected_reads=output_corrected_reads_ch
    nums_ch=bam_nums_perChr_ch
}

///////////////////////////////////////////////
/////////////COLLECTION WORKFLOWS//////////////
///////////////////////////////////////////////
workflow collect_gene_isoform_counts_perChr_wf {
  take:
    isoform_counts_ch
    gene_counts_ch
  main:

    //Collecting isoform as MTX
    isoform_mtx=collect_isoform_counts_as_mtx_perChr(isoform_counts_ch,isoform_counts_ch.map{chrom,tars,patterns,filters -> "${params.results_output}results/counts/isoform/MTX/"})
    gene_mtx=collect_gene_counts_as_mtx_perChr(gene_counts_ch,gene_counts_ch.map{chrom,tars,patterns,filters -> "${params.results_output}results/counts/gene/MTX/"})

    isoform_h5ad=collect_isoform_mtx_as_h5ad(isoform_mtx.chrom_mtx | collect, 'isoforms',"${params.results_output}results/counts/isoform/H5AD/")
    gene_h5ad=collect_gene_mtx_as_h5ad(gene_mtx.chrom_mtx | collect, 'genes',"${params.results_output}results/counts/gene/H5AD/")
  emit:
    isoform_h5ad=isoform_h5ad.h5ad_file
    isoform_mtx=isoform_mtx.chrom_mtx
    gene_h5ad=gene_h5ad.h5ad_file
    gene_mtx=gene_mtx.chrom_mtx

}
workflow collect_output_wf {

  take:
    isoform_counts_ch
    gene_counts_ch
    existing_gtf_ch
    assignment_reads_ch
    transcriptmodel_reads_ch
    corrected_reads_ch
  main:
    ///////////////////////////////////////////////////////
    //////////////////B-COLLECTING OUTPUT/////////////////////////
    ///////////////////////////////////////////////////////
    //Collect counts as MTX/H5AD
    isoform_gene_mtx_h5ad=collect_gene_isoform_counts_perChr_wf(isoform_counts_ch,gene_counts_ch)

    ///5-Collecting transcript model GTFs
    isoform_gene_mtx_h5ad.isoform_mtx.map{mtx_dir -> "${mtx_dir}/genes.tsv"}.collect().set{mtx_isoform_fs}
    existing_gtf_ch.map{chrom,gtf_tars -> gtf_tars}.collect().set{input_gtf_ch}
    gtfs=collect_gtfs(input_gtf_ch,params.gtf_f,mtx_isoform_fs,"${params.results_output}results/gtf/")
    extended_gtf=gtfs[0]
    existing_gtf=gtfs[1]
    // assignment_reads_ch.view()


}
workflow isoquant_twopass_chunked_wf {
  take:
    isoquant_preprocess_bam_perChr_ch
    chrom_sizes_f
    chunks
  main:
    ///////////////////////////////////////////////////////
    //////////////////A-FIRST PASS/////////////////////////
    ///////////////////////////////////////////////////////
    isoquant_preprocess_bam_perChr_ch.groupTuple(by:0).map{tpl -> [tpl[0],(tpl[1]).size()]}.set{bam_nums_perChr_ch}

    if (params.firstpass_mode == 'per_sample') {
        // Run isoquant once per sample on the full BAM, then fan out to all chromosomes
        per_sample_ch = isoquant_preprocess_bam_perChr_ch
            .map { chrom, sample_id, bam, bai -> tuple(sample_id, bam, bai) }
            .unique()
        per_sample_firstpass_ch = run_isoquant_firstPass_perSample(per_sample_ch, params.genedb_f, params.genome_fasta_f, params.genome_fasta_f + ".fai")
        chrom_ch_local = isoquant_preprocess_bam_perChr_ch.map { chrom, sample_id, bam, bai -> chrom }.unique()
        isoquant_firstpass_output_ch = per_sample_firstpass_ch
            .combine(chrom_ch_local)
            .map { sample_id, tar, bam, chrom -> tuple(chrom, sample_id, tar, bam) }
    } else {
        isoquant_firstpass_output_ch = run_isoquant_firstPass(isoquant_preprocess_bam_perChr_ch, params.genedb_f, params.genome_fasta_f, params.genome_fasta_f + ".fai")
    }

    // Group all samples per chromosome, then create one merged model-construction BAM per chrom
    isoquant_firstpass_output_ch
    .map{ chrom,sample_id,tar,bam -> [chrom,sample_id,tar,bam] }
    .combine(bam_nums_perChr_ch,by:0)
    .map{ chrom,sample_id,tar,bam,n -> [groupKey(chrom,n),sample_id,tar,bam] }
    .groupTuple(by:0)
    .set{ model_construction_bam_input_ch }
    model_construction_bam_ch=create_model_construction_bam_perChr(model_construction_bam_input_ch)
    // model_construction_bam_ch: [chrom, merged_bam, merged_bai]
    ////////////////////////////////////////////////////////////
    //////////////////END: A-FIRST PASS/////////////////////////
    ////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////
    //////////////////B-SECOND PASS/////////////////////////
    ///////////////////////////////////////////////////////
    ///sharding: one job per chrom on the merged BAM replaces the per-sample find + acrossSamples intersect
    mapped_unmapped_ch=find_mapped_and_unmapped_regions_perChr(model_construction_bam_ch,chrom_sizes_f)

    model_construction_bam_ch
    .combine(mapped_unmapped_ch.unmapped_bed,by:0)
    .set { modelconstructionBam_unmappedbed_tuples_ch }

    suggested_splits_ch=suggest_splits_binarySearch_merged(modelconstructionBam_unmappedbed_tuples_ch,chunks,chrom_sizes_f)


    suggested_splits_ch.flatMap { tuple ->
        def chrom = tuple[0]     // Extract chrom from the tuple
        def filePath = tuple[3]  // Get the last item in the tuple (path to file)
        file(filePath).text.split('\n') // Read file content and split into lines
            .findAll { it }             // Remove empty lines
            .collect { line ->          // Format each line and pair it with chrom
                def cols = line.split('\t') // Split line into columns by tab
                def formattedRegion = cols[0]+":"+cols[1]+"-"+cols[2] // Region formatting
                def programmaticRegion=cols[0]+"_"+cols[1]+"_"+cols[2]
                [chrom, formattedRegion,programmaticRegion] // Return chrom and formatted line as a pair
            }
    }.set {chrom_region_ch}

    /// Splitting merged BAM per chunk
    model_construction_bam_ch
      .combine(chrom_region_ch,by:0)
      .set {modelconstructionBam_Region_tuple_ch}
    modelconstructionRegionBam=split_bam_perChunk(modelconstructionBam_Region_tuple_ch,"model_construction_reads")


    ///Runnning second pass with merged single-BAM per chunk
    modelconstructionRegionBam
    .map {chrom, bam, bai, formattedRegion, programmaticRegion, counts_f  -> [chrom, bam, bai, formattedRegion, programmaticRegion ] }
    .set {isoquant_chunked_input}

    isoquant_secondpass_output_ch=run_isoquant_chunked_merged(isoquant_chunked_input, params.genedb_f, params.genome_fasta_f, params.genome_fasta_f + ".fai")

    
    isoquant_secondpass_output_ch
    ///.filter{tpl -> (tpl[1] == 'chr1_227459712_248956421') || (tpl[1] == 'chr1_159501669_170660798')}
    .set{replace_novel_names_input_ch}


    //1-Updating names of novel transcript so they don't clash between chunks
    isoquant_output_novel_names_ch=replace_novel_names(replace_novel_names_input_ch)
  

    ////////////////////////////////////////////////////////////
    //////////////////B-OUTPUT CHANNELs/////////////////////////
    ////////////////////////////////////////////////////////////
    // Each count channel emits [chrom, [tars], [patterns], [needs_chr_filter]] where patterns[i] is the
    // wildcard to extract from tars[i] and needs_chr_filter[i] signals whether rows must be filtered
    // to this chromosome (true for per-sample whole-genome first-pass tars, false for all others).
    // First-pass tars contribute N items per chrom, second-pass renamed tars contribute `chunks` — total N+chunks.
    def fp_filter = (params.firstpass_mode == 'per_sample')

    //Combining transcript-level counts channels from first and second passes
    isoquant_firstpass_output_ch
      .combine(bam_nums_perChr_ch,by:0)
      .map{chrom,sample_id,tar,bam,n -> [groupKey(chrom,n+chunks),tar,"*.transcript_grouped_tag_CB_counts.linear.tsv",fp_filter]}
      .mix(
        isoquant_output_novel_names_ch
          .combine(bam_nums_perChr_ch,by:0)
          .map{chrom,programmaticRegion,tar,n -> [groupKey(chrom,n+chunks),tar,"*.discovered_transcript_grouped_tag_CB_counts.linear.noknown.tsv",false]}
      )
      .groupTuple(by:0)
      .map{chrom,tars,patterns,filters -> [chrom,tars,patterns,filters]}
      .set{output_isoform_counts_ch}

    //Combining gene-level counts channels from first and second passes
    isoquant_firstpass_output_ch
      .combine(bam_nums_perChr_ch,by:0)
      .map{chrom,sample_id,tar,bam,n -> [groupKey(chrom,n+chunks),tar,"*.gene_grouped_tag_CB_counts.linear.tsv",fp_filter]}
      .mix(
        isoquant_output_novel_names_ch
          .combine(bam_nums_perChr_ch,by:0)
          .map{chrom,programmaticRegion,tar,n -> [groupKey(chrom,n+chunks),tar,"*.discovered_gene_grouped_tag_CB_counts.linear.tsv",false]}
      )
      .groupTuple(by:0)
      .map{chrom,tars,patterns,filters -> [chrom,tars,patterns,filters]}
      .set{output_gene_counts_ch}

    //Collecting transcript models GTFs from second-pass renamed tars
    isoquant_output_novel_names_ch
      .map{chrom,programmaticRegion,tar -> [groupKey(chrom,chunks),tar]}
      .groupTuple(by:0)
      .set{output_existing_gtf_ch}

    // corrected reads, read assignments, transcript model reads are contained in the tars
    // but not consumed by any downstream collection process, so emit empty channels
    output_corrected_reads_ch = Channel.empty()
    output_assignment_reads_ch = Channel.empty()
    output_transcriptmodel_reads_ch = Channel.empty()


    ///////////////////////////////////////////////////////
    //////////////////END: OUTPUT CHANNELS/////////////////
    ///////////////////////////////////////////////////////

  emit:
    isoform_counts=output_isoform_counts_ch
    gene_counts=output_gene_counts_ch
    existing_gtf=output_existing_gtf_ch
    assignment_reads=output_assignment_reads_ch
    transcriptmodel_reads=output_transcriptmodel_reads_ch
    corrected_reads=output_corrected_reads_ch
    nums_ch=bam_nums_perChr_ch

}


workflow ISOQUANT_TWOPASS_WF {

  take:
    fullBam_ch
  main:
    if (fullBam_ch == 'independent workflow'){
      def samplesheet_header = file(params.input_samples_path).readLines()[0].split(',') as List
      if ('mapped_bam' in samplesheet_header) {
        Channel
          .fromPath(params.input_samples_path)
          .splitCsv(sep: ',', header: true)
          .map { it -> [it.sample_id, it.mapped_bam, it.mapped_bam + ".bai"] }
          .set { fullBam_ch }
      } else {
        Channel
          .fromPath(params.input_samples_path)
          .splitCsv(sep: ',', header: true)
          .map { it ->
            def sample_id = it.sample_id
            def bam_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam"
            def bai_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam.bai"
            [sample_id, bam_path, bai_path]
          }
          .set { fullBam_ch }
        Channel
          .fromPath(params.input_samples_path)
          .splitCsv(sep: ',', header: true)
          .map { it -> [it.sample_id, it.nr_samples_multiplexed] }
          .set { sampleNames_nrDons }

        fullBam_ch =  fullBam_ch.combine(sampleNames_nrDons, by: 0)
        fullBam_ch=fullBam_ch
        .filter { experiment, bam, bai, npooled -> npooled == '1' }
        .map { experiment, bam, bai, npooled -> [experiment, bam, bai] }

        def deconv_dir_path = file("${params.results_output}deconvolution/bam/")
        if (deconv_dir_path.exists()) {
          Channel
            .fromPath(params.input_samples_path)
            .splitCsv(sep: ',', header: true)
            .filter { it -> it.nr_samples_multiplexed != '1' }
            .flatMap { it ->
              def sample_id = it.sample_id
              def nr_donors = it.nr_samples_multiplexed as int
              (0..<nr_donors).collect { n ->
                def bam_path = "${params.results_output}deconvolution/bam/${sample_id}_donor${n}.nosupplementary.bam"
                tuple("${sample_id}_donor${n}" as String, file(bam_path), file("${bam_path}.bai"))
              }
            }
            .set { split_bams }
          fullBam_ch = fullBam_ch.mix(split_bams)
        }
      }
    }
    def chromosomes_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                            'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                            'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY','chrM']

    if(params.chunks) {
    chrom_sizes_f=params.nf_basedir+"data/hg38.chrom.sizes"
    }

    chrom_ch=chroms(chromosomes_list)

    preprocessed_bam_perChr_ch=chrom_ch.combine(fullBam_ch)

    //Processing chrM separately
    preprocessed_bam_perChr_ch.filter{chrom,sample_id,bam,bai -> chrom=="chrM"}.set{preprocessed_bam_chrM_ch}
    chrM_output_chs=isoquant_chrM(preprocessed_bam_chrM_ch,chrom_sizes_f)

    //Processing other chromosomes separately
    preprocessed_bam_perChr_ch.filter{chrom,sample_id,bam,bai -> chrom!="chrM"}.set{preprocessed_bam_nochrM_ch}
    nochrM_output_chs=isoquant_twopass_chunked_wf(preprocessed_bam_nochrM_ch,chrom_sizes_f,params.chunks)

    //Collecting output: MTX, GTF, reads
    collect_output_wf(
    (nochrM_output_chs.isoform_counts).concat(chrM_output_chs.isoform_counts),
    (nochrM_output_chs.gene_counts).concat(chrM_output_chs.gene_counts),
    (nochrM_output_chs.existing_gtf).concat(chrM_output_chs.existing_gtf),
    (nochrM_output_chs.assignment_reads).concat(chrM_output_chs.assignment_reads),
    (nochrM_output_chs.transcriptmodel_reads).concat(chrM_output_chs.transcriptmodel_reads),
    (nochrM_output_chs.corrected_reads).concat(chrM_output_chs.corrected_reads)
    )
}
