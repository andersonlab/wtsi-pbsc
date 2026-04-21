nextflow.enable.dsl=2

///Modules
include {SQANTI3_QC; SQANTI3_FILTER} from './modules/sqanti3.nf'

///Subworkflows
include {BAM_PROCESSING} from './subworkflows/bam_processing/bam_processing.nf'
include {DECONVOLUTION} from './subworkflows/deconvolution/deconvolution.nf'
include {ISOQUANT_TWOPASS_PROCESS} from './subworkflows/isoquant_recipes/isoquant_twopass_process.nf'

include {mtx_subset_wf} from './subworkflows/core/mtx_subset.nf'


//include {customPublish as customPublishFilteredH5ADIsoform} from './modules/customPublish.nf'
//include {customPublish as customPublishFilteredMTXIsoform} from  './modules/customPublish.nf'


/// Setting default parameters
if(!params.barcode_correction_percentile) {
  params.barcode_correction_percentile=98
}
if(!params.min_polya_length) {
  params.min_polya_length=20
}

//if(!params.exclude_samples) {
//  params.exclude_samples=[]
//}


assert params.run_mode in ['with_quant', 'pre_quant'] : "ERROR: params.run_mode must be one of: 'with_quant', 'pre_quant'"


workflow full{
  BAM_PROCESSING()
  mapped_reads=BAM_PROCESSING.out.mapped_reads
  if (params.run_deconvolution == 'TRUE') {
    DECONVOLUTION(mapped_reads)
    mapped_reads=DECONVOLUTION.out.fullBam_ch
  }
  mapped_reads.view()
  if (params.run_mode == 'with_quant'){ 
    ISOQUANT_TWOPASS_PROCESS(mapped_reads)
  }
}

workflow bam_processing_wf {
    // Independent workflow entry for deconvolution
    BAM_PROCESSING()
}

workflow deconvolution_wf {
    // Independent workflow entry for deconvolution
    input_ch = 'independent workflow'
    DECONVOLUTION(input_ch)
}

workflow isoquant_twopass_wf {
    input_ch = 'independent workflow'
    ISOQUANT_TWOPASS_PROCESS(input_ch)
}


///////////////////////////////////////
//recheck below
///////////////////////////////////////



workflow sqanti3 {

  def input_gtf_f="${params.results_output}results/gtf/transcript_models.gtf"
  sqanti3_qc_ch=SQANTI3_QC(input_gtf_f,params.gtf_f,params.genome_fasta_f,params.polya_f,params.cage_peak_f,params.polya_sites,"/nfs/team152/oe2/isogut/software/SQANTI3-5.3.0/")
  sqanti3_filter_ch=SQANTI3_FILTER(sqanti3_qc_ch.map{output_dir -> "${output_dir}/transcript_models_classification.txt"},input_gtf_f,params.sqanti3_path)

  sqanti3_qc_ch.view()
  Channel
  .fromPath("${params.results_output}results/counts/isoform/MTX/*/matrix.mtx")
  .map{path -> path.parent}
  .set{prefiltered_mtx_dir_ch}
  mtx_subset_output_ch=mtx_subset_wf(prefiltered_mtx_dir_ch,sqanti3_filter_ch.map{sqanti3_filter_dir -> "${sqanti3_filter_dir}/transcript_models_inclusion-list.txt"})
  mtx_subset_output_ch.h5ad_file.view()

///  customPublishFilteredH5ADIsoform(mtx_subset_output_ch.h5ad_file,"${params.results_output}results/counts_sqanti3/isoform/H5AD/")
///  customPublishFilteredMTXIsoform((mtx_subset_output_ch.isoform_mtx).collect(),"${params.results_output}results/counts_sqanti3/isoform/MTX/")

}

