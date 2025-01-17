nextflow.enable.dsl=2


include { split_reads; remove_primer; tag_bam; refine_reads } from './modules/fltnc.nf'
include { barcode_correction; dedup_reads } from './modules/barcodes.nf'
include { pbmm2 } from './modules/pbmm2.nf'
include { isoquant_split_by_chr; run_isoquant } from './modules/isoquant.nf'


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
  dedup               = dedup_reads(barcode_corrected.barcode_corrected_tuple,params.barcode_correction_method,params.barcode_correction_percentile)

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
    def bam_path = "${params.results_output}qc/mapped/${sample_id}.dedup.bam"
    [sample_id, bam_path]
  }
  .set { dedup_bam_tuples }
  dedup_bam_tuples.view()
  mapped_reads         = pbmm2(dedup_bam_tuples,params.genome_fasta_f)
}

workflow isoquant {


def chromosomes_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                        'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

def excluded_samples_list = file(params.exclude_samples).exists() ?
  file(params.exclude_samples).text.split('\n').drop(1).collect { it.split(',')[0].trim() } : []


Channel
  .fromPath(params.input_samples_path)
  .splitCsv(sep: ',', header: true)
  .filter { it -> !(it.sample_id in excluded_samples_list) }
  .flatMap { it ->
    def sample_id = it.sample_id
    def bam_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam"
    def bai_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam.bai"
    chromosomes_list.collect { chromosome -> [sample_id, chromosome, bam_path, bai_path] }
  }
  .set { mapped_bam_chrom_tuples }


isoquant_split_tuple=isoquant_split_by_chr(mapped_bam_chrom_tuples)

isoquant_split_tuple
    .map { sample_id, chrom, bam, bai -> [chrom, sample_id, bam, bai] }
    .groupTuple()
    .set { bam_chrom_tuples }

isoquant_output=run_isoquant(bam_chrom_tuples,params.gtf_f,params.genome_fasta_f)

}

workflow collapse_classify {

    def mapped_reads_paths = "${params.results_output}/results/mapped_bam/*mapped.realcells_only.bam"
    Channel.fromPath(mapped_reads_paths).set { mapped_bams }
    collapsed_isoforms   = collapse_isoforms(mapped_bams)

    //Preparing input for SQANTI3 QC
    collapsed_isoforms_tuple = collapsed_isoforms.collapsed_gff_path.map { path ->
       def sampleId = path.name.find(/(Isogut\d+)/)
       tuple(sampleId, path)
    }
    collapsed_abundance_tuple = collapsed_isoforms.collapsed_abundance_path.map { path ->
       def sampleId = path.name.find(/(Isogut\d+)/)
       tuple(sampleId, path)
    }

    collapsed_isoforms_tuple
       .join(collapsed_abundance_tuple, by: 0)
       .set { sqanti3_qc_input }
    sqanti3_qc = sqanti3_classify_isoforms(sqanti3_qc_input)

    //Preparing input for SQANTI3 FILTERING
    sqanti3_classification_tuple = sqanti3_qc.classification_files.map { path ->
       def sampleId = path.name.find(/(Isogut\d+)/)
       tuple(sampleId, path)
    }
    sqanti3_faa_tuple = sqanti3_qc.faa_files.map { path ->
       def sampleId = path.name.find(/(Isogut\d+)/)
       tuple(sampleId, path)
    }

    sqanti3_gtf_tuple = sqanti3_qc.gtf_files.map { path ->
       def sampleId = path.name.find(/(Isogut\d+)/)
       tuple(sampleId, path)
    }

    sqanti3_classification_tuple
       .join(sqanti3_gtf_tuple, by: 0)
       .join(sqanti3_faa_tuple, by: 0)
       .set { sqanti3_filter_input }
     sqanti3_filter = sqanti3_filter_isoforms(sqanti3_filter_input)

     ///////////

    def fasta_paths="${params.results_output}/results/unmapped_bam/*.fasta"
    fastaFiles = Channel
       .fromPath(fasta_paths)
       .map { file ->
           // Extract just the sample identifier from the filename
           def sampleId = file.name.find(/(Isogut\d+)/)
           tuple(sampleId, file)
    }
    sqanti3_classification_tuple = sqanti3_qc.classification_files.map { path ->
       def sampleId = path.name.find(/(Isogut\d+)/)
       tuple(sampleId, path)
    }
    sqanti3_filtered_tuple = sqanti3_filter.filtered_classification_files.map { path ->
       def sampleId = path.name.find(/(Isogut\d+)/)
       tuple(sampleId, path)
    }

   collapsed_group_path = collapsed_isoforms.collapsed_group_path.map { path ->
       def sampleId = path.name.find(/(Isogut\d+)/)
       tuple(sampleId, path)
    }
    fastaFiles
        .join(sqanti3_classification_tuple, by: 0)
        .join(sqanti3_filtered_tuple, by: 0)
        .join(collapsed_abundance_tuple, by: 0)
        .join(collapsed_group_path, by: 0)
        .set { sqanti3_seurat_input }
    sqanti3_seurat = sqanti3_seurat_isoform(sqanti3_seurat_input)


    //prepared_gff         = prepare_input_GFF(collapsed_isoforms.collapsed_gff_path)
    //classified_isoforms  = classify_isoforms(prepared_gff,collapsed_isoforms.collapsed_abundance_path)
    // filtered_isoforms    = filter_isoforms(classified_isoforms.classification_files,classified_isoforms.junction_files,collapsed_isoforms.collapsed_gff_path)

    // Preparing all required files for seurat_isoform
    // Extract the identifier from the fasta filenames
    // def fasta_paths="${params.results_output}/results/unmapped_bam/*.fasta"
    // fastaFiles = Channel
    //    .fromPath(fasta_paths)
    //    .map { file ->
    //        // Extract just the sample identifier from the filename
    //        def sampleId = file.name.find(/(Isogut\d+)/)
    //        tuple(sampleId, file)
    //    }

    // Setup for the classificationFiles channel
    // classificationFiles = filtered_isoforms.classification_filtered_files.map { path ->
    //    def sampleId = path.name.find(/(Isogut\d+)/)
    //    tuple(sampleId, path)
    // }

    // Setup for the collapsed_group_path channel
    // collapsed_group_path = collapsed_isoforms.collapsed_group_path.map { path ->
    //    def sampleId = path.name.find(/(Isogut\d+)/)
    //    tuple(sampleId, path)
    // }

    // Join channels by matching the sample identifier, which is the base name transformation
    // matchedFiles = fastaFiles
    //    .join(classificationFiles, by: 0)
    //    .join(collapsed_group_path, by: 0)
    //    .set { seurat_input }

    // seurat_output        = seurat_isoform(seurat_input)


}
