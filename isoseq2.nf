nextflow.enable.dsl=2


///Modules
include { split_reads; remove_primer; tag_bam; refine_reads } from './modules/fltnc.nf'
include {barcode_correction; get_barcodes; supset_bam; dedup_reads; combine_dedups; bam_stats} from './modules/barcodes.nf'

include { pbmm2 } from './modules/pbmm2.nf'
include {SQANTI3_QC; SQANTI3_FILTER} from './modules/sqanti3.nf'
///include {preprocess_bam; find_mapped_and_unmapped_regions_per_sampleChrom; acrossSamples_mapped_unmapped_regions_perChr; suggest_splits_binarySearch; split_bams; create_genedb_fasta_perChr; run_isoquant_chunked} from './modules/isoquant.nf'


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
    /// refined_bam_stats = postrefine_stats(refined_reads.refined_bam)
}



workflow correct_barcodes {

def excluded_samples_list = file(params.exclude_samples).exists() ?
  file(params.exclude_samples).text.split('\n').drop(1).collect { it.split(',')[0].trim() } : []

println "Excluded samples: $excluded_samples_list"

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

barcode_corrected   = barcode_correction(refined_bam_tuples,params.threeprime_whitelist,params.barcode_correction_method,params.barcode_correction_percentile)

get_barcodes(barcode_corrected.barcode_corrected_tuple, params.number_of_chunks)
  ////
barcode_channel=get_barcodes.out.barcodes_tuple.transpose()

combined_ch = barcode_corrected.barcode_corrected_tuple.combine(barcode_channel, by: 0)
//combined_ch = barcode_corrected.combine(barcode_channel, by: 0)

supset_bam(combined_ch)

dedup_reads(supset_bam.out.chunk_tuple, params.dedup_batch_size)

deduped_chunks_ch=dedup_reads.out.dedup_tuple.groupTuple()

combine_dedups(deduped_chunks_ch)

combined_dedup_ch=combine_dedups.out.dedup_tuple

bam_stats(combined_dedup_ch, params.barcode_correction_method,params.barcode_correction_percentile)

}


workflow map_pbmm2 {

def excluded_samples_list = file(params.exclude_samples).exists() ?
  file(params.exclude_samples).text.split('\n').drop(1).collect { it.split(',')[0].trim() } : []

println "Excluded samples: $excluded_samples_list"

Channel
  .fromPath(params.input_samples_path)
  .splitCsv(sep: ',', header: true)
  .filter { it -> !(it.sample_id in excluded_samples_list) }
  .map { it ->
    def sample_id = it.sample_id
    def bam_path = "${params.results_output}qc/dedup/${sample_id}.dedup.bam"//mapped to dedup
    [sample_id, bam_path]
  }
  .set { dedup_bam_tuples }
  dedup_bam_tuples.view()
  mapped_reads         = pbmm2(dedup_bam_tuples,params.genome_fasta_f)
}


workflow isoquant_twopass {


  def chromosomes_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                          'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                          'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY','chrM']
  ///def chromosomes_list = ['chr22']

  chrom_ch=chroms(chromosomes_list)
  ///.filter{chrom -> chrom=='chr2' }
  fullBam_ch=bamsWithExclusion()
  ///.filter{tpl -> (tpl[0]=='Isogut14548280') || (tpl[0]=='Isogut14548279') || (tpl[0]=='Isogut14548278') || (tpl[0]=='Isogut14548277') || (tpl[0]=='Isogut14548276') || (tpl[0]=='Isogut14548275')}
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

workflow isoquant_chunked_old {


def chromosomes_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                        'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

///Testing
/// def chromosomes_list = ['chr22']
///
chrom_ch=chroms(chromosomes_list)
fullBam_ch=bamsWithExclusion()

///Similar to fullBam_ch, but emits processed BAMs
processedFullBam_ch=preprocess_bam(fullBam_ch)
///Similar to processedFullBam_ch, but has a record per chromosome (useful for processing by chromosome using the full processed BAM files)
///TESTING: .filter{record -> (record[0]=='Isogut15020393') || (record[0]=='Isogut14548284') }
processedFullBam_ch
.combine(chrom_ch)
.map{ sample_id,bam,bai,chrom -> [sample_id,chrom,bam,bai] }
.set{processedFullBam_Chrom_tuple_ch}

///Similar to processedFullBam_Chrom_tuple_ch but tuples are grouped per chromosome
processedFullBam_Chrom_tuple_ch
.map{sample_id,chrom,bam,bai -> [chrom,sample_id,bam,bai]}
.groupTuple()
.set{processedFullBam_Chrom_groupedTuple_ch}

///////////////////////////////////
///////////TO BE REMOVED///////////
///////////////////////////////////
///mapped_bam_ch
/// .flatMap { it ->
///      def sample_id = it.sample_id
///      def bam_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam"
///      def bai_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam.bai"
///      chromosomes_list.collect { chromosome -> [sample_id, chromosome, bam_path, bai_path] }
///    }
///    .set { mapped_bam_chrom_tuples }
///////////////////////////////////
///////////END: TO BE REMOVED//////
///////////////////////////////////







    /////////////////////////////////////////////////////////
    ///////////Workflow for splitting by chromosome only/////
    /////////////////////////////////////////////////////////


    ////// isoquant_split_tuple=isoquant_split_by_chr(mapped_bam_chrom_tuples)

    ////// isoquant_split_tuple
    //////    .map { sample_id, chrom, bam, bai -> [chrom, sample_id, bam, bai] }
    //////    .groupTuple()
    //////    .set { bam_chrom_tuples }

    ////// isoquant_output=run_isoquant(bam_chrom_tuples,params.gtf_f,params.genome_fasta_f)

    //////////////////////////////////////////////////////////////////////
    /////////Workflow for splitting by chromosome only and chunks/////////
    //////////////////////////////////////////////////////////////////////

    def unmapped_regions_f="${params.results_output}misc/unmapped_regions.bed"
    def mapped_regions_f="${params.results_output}misc/mapped_regions.bed"
    def chunks_numreads_f="${params.results_output}misc/chunks_numreads.csv"




    ///Finding mapped/unmapped regions (per chromosome)
    mapped_unmapped_regions_tuple_ch=find_mapped_and_unmapped_regions_per_sampleChrom(processedFullBam_Chrom_tuple_ch,chrom_sizes_f)
    mapped_unmapped_regions_tuple_ch.groupTuple().set{mapped_unmapped_regions_groupedTuple_ch}

    ///Finding intersections of unmapped regions, and merging mapped regions (per chromosome)
    acrossSamples_mapped_unmapped_regions_perChr_ch=acrossSamples_mapped_unmapped_regions_perChr(mapped_unmapped_regions_groupedTuple_ch)

    /// Merging grouped BAM tuples with mapped/unmapped BED files
    processedFullBam_Chrom_groupedTuple_ch.join(acrossSamples_mapped_unmapped_regions_perChr_ch.mapped_bed,by:0).set { processedFullBam_mappedbed_tuples_ch }
    processedFullBam_Chrom_groupedTuple_ch.join(acrossSamples_mapped_unmapped_regions_perChr_ch.unmapped_bed,by:0).set { processedFullBam_unmappedbed_tuples_ch }
    suggested_splits_ch=suggest_splits_binarySearch(processedFullBam_unmappedbed_tuples_ch,params.chunks,chrom_sizes_f)

    ///Logging mapped/unmapped BEDs
    ///acrossSamples_mapped_unmapped_regions_perChr.mapped_bed
    ///    .collectFile(name: mapped_regions_f, cache: false)
    //////    .view(f -> f.toString())

    ///acrossSamples_mapped_unmapped_regions_perChr.unmapped_bed
    ///    .collectFile(name: unmapped_regions_f, cache: false)
    //////    .view(f -> f.toString())



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

    /// Splitting BAMs according to suggested regions
    processedFullBam_Chrom_groupedTuple_ch
      .combine(chrom_region_ch,by:0)
      .set {processedFullBam_Region_groupedTuple_ch}
    processedRegionBam=split_bams(processedFullBam_Region_groupedTuple_ch)
    /// Logging numreads per chunk to be analyzed prior to running chunks (to better understand memory allocation)
    processedRegionBam
    .map {chrom, sample_ids, bams, bais, formattedRegion, programmaticRegion, counts  -> counts }
    .collectFile(name: chunks_numreads_f, cache: false)
    //////.view(f -> f.toString())



    /// Running chunked isoquant

    /// Splitting Genome FASTA and GENCODE GTF per chromosome
    chrom_ch
    .map {chrom -> [chrom,params.gtf_f,params.genome_fasta_f]}
    .set {chrom_genedb_fasta_chr_input_ch}
    chrom_genedb_fasta_chr_ch=create_genedb_fasta_perChr(chrom_genedb_fasta_chr_input_ch)


    processedRegionBam
    .map {chrom, sample_ids, bams, bais, formattedRegion, programmaticRegion, counts  -> [chrom, sample_ids, bams, bais, formattedRegion, programmaticRegion, file(counts).splitCsv(header:false)[0][1] ] }
    .combine(chrom_genedb_fasta_chr_ch.map{chrom, chrom_gtf_f, chrom_genedb_f,chrom_fasta, chrom_fai -> [chrom,chrom_genedb_f,chrom_fasta,chrom_fai]},by:0)
    .set {isoquant_chunked_input}


    isoquant_output=run_isoquant_chunked(isoquant_chunked_input)
    isoquant_output.view()

}
