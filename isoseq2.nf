nextflow.enable.dsl=2

///Modules
include { split_reads; remove_primer; tag_bam; refine_reads } from './modules/fltnc.nf'
include {barcode_correction; get_barcodes; supset_bam; supset_bam_with_bai; dedup_reads; combine_dedups; combine_mupped; bam_stats} from './modules/barcodes.nf'

include { pbmm2 } from './modules/pbmm2.nf'
include {SQANTI3_QC; SQANTI3_FILTER} from './modules/sqanti3.nf'
///include {preprocess_bam; find_mapped_and_unmapped_regions_per_sampleChrom; acrossSamples_mapped_unmapped_regions_perChr; suggest_splits_binarySearch; split_bams; create_genedb_fasta_perChr; run_isoquant_chunked} from './modules/isoquant.nf'
include { mpileup; cellsnp; vireo; subset_vcf }  from './modules/deconvolution.nf'

///Subworkflows
include {chroms} from './subworkflows/core/chroms.nf'
include {bamsWithExclusion} from './subworkflows/core/bamsWithExclusion.nf'
include {preprocess_bam_perChr_wf} from './subworkflows/isoquant_components/preprocess_bam.nf'
include {genedb_perChr_wf} from './subworkflows/isoquant_components/genedb.nf'
include {isoquant_twopass_perChr_wf; isoquant_twopass_chunked_wf; isoquant_chrM;collect_output_wf} from './subworkflows/isoquant_recipes/isoquant_twopass.nf'
include {mtx_subset_wf} from './subworkflows/core/mtx_subset.nf'


include {customPublish as customPublishFilteredH5ADIsoform} from './modules/customPublish.nf'
include {customPublish as customPublishFilteredMTXIsoform} from  './modules/customPublish.nf'


/// Setting default parameters
if(!params.barcode_correction_percentile) {
  params.barcode_correction_percentile=99
}
if(!params.min_polya_length) {
  params.min_polya_length=20
}

if(!params.exclude_samples) {
  params.exclude_samples=[]
}

if(params.chunks) {
  chrom_sizes_f=params.nf_basedir+"data/hg38.chrom.sizes"
}

/// This QC workflow will do the pre-processing in: https://isoseq.how/umi/cli-workflow.html
workflow fltnc {
    main:
      /// Obtaining (sample_id,bam_file) tuples from the input_samples.csv file
      Channel
        .fromPath(params.input_samples_path)
        .splitCsv(sep: ',', header: true)
        .map { it -> [it.sample_id, it.long_read_path] } // Create a tuple with bam_path and sample_id
        .set { hifi_bam_tuples }

<<<<<<< HEAD
      /// Every process from now on outputs a (sample_id,bam_file) tuple which is fed on to the next process
      reads_split           = split_reads(hifi_bam_tuples,params.skera_primers)
      primer_removed        = remove_primer(reads_split.split_reads_tuple,params.tenx_primers)
      tagged                = tag_bam(primer_removed.removed_primer_tuple)
      refined_bam_tuples         = refine_reads(tagged.tagged_tuple,params.tenx_primers,params.min_polya_length)
      /// refined_bam_stats = postrefine_stats(refined_reads.refined_bam)
    emit:
      refined_bam_tuples
=======
    /// Obtaining (sample_id,bam_file) tuples from the input_samples.csv file
    Channel
    	.fromPath(params.input_samples_path)
    	.splitCsv(sep: ',', header: true)
    	.map { it -> [it.sample_id, it.long_read_path] } // Create a tuple with bam_path and sample_id
    	.set { hifi_bam_tuples }

    /// Every process from now on outputs a (sample_id,bam_file) tuple which is fed on to the next process
    reads_split           = split_reads(hifi_bam_tuples,params.skera_primers)
    primer_removed        = remove_primer(reads_split.split_reads_tuple,params.tenx_primers)
    tagged                = tag_bam(primer_removed.removed_primer_tuple)
    refined_reads         = refine_reads(tagged.tagged_tuple,params.tenx_primers,params.min_polya_length)
    //refined_bam_stats = postrefine_stats(refined_reads.refined_bam)
>>>>>>> origin/main
}



workflow correct_barcodes_process {
  take:
    refined_bam_tuples
  main:
      def excluded_samples_list = file(params.exclude_samples).exists() ?
        file(params.exclude_samples).text.split('\n').drop(1).collect { it.split(',')[0].trim() } : []

      println "Excluded samples: $excluded_samples_list"
      if (refined_bam_tuples == 'independent workflow'){
          Channel
            .fromPath(params.input_samples_path)
            .splitCsv(sep: ',', header: true)
            .filter { it -> !(it.sample_id in excluded_samples_list) }
            .map { it ->
              def sample_id = it.sample_id
              def bam_path = "${params.results_output}qc/refined/${sample_id}.fltnc.bam"
              [sample_id, bam_path]
            }
            .set { refined_bam_tuples }
      }

      barcode_corrected   = barcode_correction(refined_bam_tuples,params.threeprime_whitelist,params.barcode_correction_method,params.barcode_correction_percentile)

      get_barcodes(barcode_corrected.barcode_corrected_tuple, params.number_of_chunks)
      
      barcode_channel=get_barcodes.out.barcodes_tuple.transpose()

      combined_ch = barcode_corrected.barcode_corrected_tuple.combine(barcode_channel, by: 0)

      supset_bam(combined_ch)

      dedup_reads(supset_bam.out.chunk_tuple, params.dedup_batch_size)

      mapped_chunks = pbmm2(dedup_reads.out.dedup_tuple, params.genome_fasta_f)

      mapped_chunks_ch=mapped_chunks.map_tuple.groupTuple()

      mapped_reads=combine_mupped(mapped_chunks_ch)

      //combined_dedup_ch=combine_dedups.out.dedup_tuple

      //dedup_bam_tuples = bam_stats(combined_dedup_ch, params.barcode_correction_method,params.barcode_correction_percentile)

  emit:
    mapped_reads

}


workflow correct_barcodes {
    // Independent workflow entry for isoquant_twopass
    input_ch = 'independent workflow'
    correct_barcodes_process(input_ch)
}

workflow isoquant_twopass_process {

  take:
    fullBam_ch
  main:

    def chromosomes_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                            'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                            'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY','chrM']

    chrom_ch=chroms(chromosomes_list)

    if (fullBam_ch == 'independent workflow'){
      fullBam_ch=bamsWithExclusion()
    }

    chrom_genedb_fasta_chr_ch=genedb_perChr_wf(chrom_ch,params.gtf_f,params.genome_fasta_f)
    preprocessed_bam_perChr_ch=preprocess_bam_perChr_wf(chrom_ch,fullBam_ch)

    //Processing chrM separately
    preprocessed_bam_perChr_ch.filter{chrom,sample_id,bam,bai -> chrom=="chrM"}.set{preprocessed_bam_chrM_ch}
    chrom_genedb_fasta_chr_ch.filter{chrom,gene_db,fa,fai -> chrom=="chrM"}.set{chrom_genedb_fasta_chrM_ch}
    chrM_output_chs=isoquant_chrM(preprocessed_bam_chrM_ch,chrom_genedb_fasta_chrM_ch,chrom_sizes_f)

    //Processing other chromosomes separately
    preprocessed_bam_perChr_ch.filter{chrom,sample_id,bam,bai -> chrom!="chrM"}.set{preprocessed_bam_nochrM_ch}
    chrom_genedb_fasta_chr_ch.filter{chrom,gene_db,fa,fai -> chrom!="chrM"}.set{chrom_genedb_fasta_nochrM_ch}
    nochrM_output_chs=isoquant_twopass_chunked_wf(preprocessed_bam_nochrM_ch,chrom_genedb_fasta_nochrM_ch,chrom_sizes_f,params.chunks)

    //Collecting output: MTX, GTF, reads
    collect_output_wf(
    (nochrM_output_chs.isoform_counts).concat(chrM_output_chs.isoform_counts),
    (nochrM_output_chs.gene_counts).concat(chrM_output_chs.gene_counts),
    (nochrM_output_chs.existing_gtf).concat(chrM_output_chs.existing_gtf),
    (nochrM_output_chs.extended_gtf).concat(chrM_output_chs.extended_gtf),
    (nochrM_output_chs.assignment_reads).concat(chrM_output_chs.assignment_reads),
    (nochrM_output_chs.transcriptmodel_reads).concat(chrM_output_chs.transcriptmodel_reads),
    (nochrM_output_chs.corrected_reads).concat(chrM_output_chs.corrected_reads)
    )
}

  chrom_ch=chroms(chromosomes_list)
  ///.filter{chrom -> chrom=='chr2' }
  fullBam_ch=bamsWithExclusion()
  ///.filter{tpl -> (tpl[0]=='Isogut14548280') || (tpl[0]=='Isogut14548279') || (tpl[0]=='Isogut14548278') || (tpl[0]=='Isogut14548277') || (tpl[0]=='Isogut14548276') || (tpl[0]=='Isogut14548275')}
  chrom_genedb_fasta_chr_ch=genedb_perChr_wf(chrom_ch,params.gtf_f,params.genome_fasta_f)
  preprocessed_bam_perChr_ch=preprocess_bam_perChr_wf(chrom_ch,fullBam_ch)

  //TODO: implement isoquant_twopass_perChr_wf
}
workflow isoquant_twopass {
    // Independent workflow entry for isoquant_twopass
    input_ch = 'independent workflow'
    isoquant_twopass_process(input_ch)
}

workflow sqanti3 {


  def input_gtf_f="${params.results_output}results/gtf/transcript_models.gtf"
  sqanti3_qc_ch=SQANTI3_QC(input_gtf_f,params.gtf_f,params.genome_fasta_f,params.polya_f,params.cage_peak_f,params.polya_sites,"/nfs/team152/oe2/isogut/software/SQANTI3-5.3.0/")
  sqanti3_filter_ch=SQANTI3_FILTER(sqanti3_qc_ch.map{output_dir -> "${output_dir}/transcript_models_classification.txt"},input_gtf_f,params.sqanti3_path)

  sqanti3_qc_ch.view()
  Channel
  .fromPath("${params.results_output}results/counts/isoform/MTX/*/matrix.mtx")
  .map{it -> it.parent}
  .set{prefiltered_mtx_dir_ch}
  mtx_subset_output_ch=mtx_subset_wf(prefiltered_mtx_dir_ch,sqanti3_filter_ch.map{sqanti3_filter_dir -> "${sqanti3_filter_dir}/transcript_models_inclusion-list.txt"})
  mtx_subset_output_ch.h5ad_file.view()

///  customPublishFilteredH5ADIsoform(mtx_subset_output_ch.h5ad_file,"${params.results_output}results/counts_sqanti3/isoform/H5AD/")
///  customPublishFilteredMTXIsoform((mtx_subset_output_ch.isoform_mtx).collect(),"${params.results_output}results/counts_sqanti3/isoform/MTX/")

}

workflow deconvolution{
  take:
    mapped_reads
  main:

    def excluded_samples_list = file(params.exclude_samples).exists() ?
      file(params.exclude_samples).text.split('\n').drop(1).collect { it.split(',')[0].trim() } : []

    if (mapped_reads == 'independent workflow'){
    Channel
      .fromPath(params.input_samples_path)
      .splitCsv(sep: ',', header: true)
      .filter { it -> !(it.sample_id in excluded_samples_list) }
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

    mpileup_out_chanel = mpileup(for_deconv,params.genome_fasta_f)
    mpileup_out_chanel_subset = subset_vcf(mpileup_out_chanel,params.subset_regions_bed)
    cellsnp_out = cellsnp(mpileup_out_chanel_subset)

    sampleNames_cellsnp_nrDons = cellsnp_out.combine(sampleNames_nrDons, by: 0)
    vireo(sampleNames_cellsnp_nrDons)
    
    barcode_channel=vireo.out.barcodes_tuple.transpose()
    combined_ch = for_deconv.combine(barcode_channel, by: 0)
    supset_bam_with_bai(combined_ch)
    fullBam_ch_pre = supset_bam_with_bai.out.per_donor_tuple

    fullBam_ch_pre.mix(not_for_deconv).set{fullBam_ch}
    
  emit:
    fullBam_ch

}

workflow deconvolutionwf {
    // Independent workflow entry for deconvolution
    input_ch = 'independent workflow'
    deconvolution(input_ch)
}

workflow full{

  fltnc()
  correct_barcodes_process(fltnc.out.refined_bam_tuples)
  deconvolution(correct_barcodes_process.out.mapped_reads)
  isoquant_twopass_process(deconvolution.out.fullBam_ch)
}
