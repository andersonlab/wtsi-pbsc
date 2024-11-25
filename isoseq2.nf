nextflow.enable.dsl=2


include { split_reads; remove_primer; tag_bam; refine_reads } from './modules/fltnc.nf'
/// include { barcodes } from './modules/barcodes.nf'
/// include { dedup } from './modules/dedup.nf'
/// include { pbmm2 } from './modules/pbmm2.nf'







/// This QC workflow will do the pre-processing in: https://isoseq.how/umi/cli-workflow.html
workflow fltnc {

    /// Obtaining (sample_id,bam_file) tuples from the input_samples.csv file
    Channel
    	.fromPath(params.input_samples_path)
    	.splitCsv(sep: ',', header: true)
    	.map { it -> [it.sample_id, it.long_read_path] } // Create a tuple with bam_path and sample_id
    	.set { sample_id_and_bam }

    /// Every process from now on outputs a (sample_id,bam_file) tuple which is fed on to the next process
    reads_split           = split_reads(sample_id_and_bam,params.skera_primers)
    primer_removed        = remove_primer(reads_split.split_reads_tuple,params.tenx_primers)
    tagged                = tag_bam(primer_removed.removed_primer_tuple)
    refined_reads         = refine_reads(tagged.tagged_tuple,params.tenx_primers,params.min_polya_length)
    /// refined_bam_stats = postrefine_stats(refined_reads.refined_bam)
}



workflow correct_dedup {

 def refined_reads_paths = "${params.results_output}/qc/refined/*.bam"

  Channel
    .fromPath(refined_reads_paths)
    .set {refined_bams}

  sorted_barcodes   = barcode_correction(refined_bams)
  sorted_cells      = sort_cellbarcodes(sorted_barcodes.barcode_corrected_bam)
  deduplication     = dedup_reads(sorted_cells.sorted_bam)
  stats             = postdedup_stats(deduplication.dedup_bam)
}

workflow map_to_genome {

    def unmapped_reads_paths = "${params.results_output}/results/unmapped_bam/*.bam"
    Channel.fromPath(unmapped_reads_paths).set { unmapped_bams }
    mapped_reads         = map_reads(unmapped_bams)
    mapped_realcells_only = make_real_cells_bams(mapped_reads.mapped_bam)



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

workflow isoquant {
    // Create a channel to collect all BAM files
    Channel
        .fromPath("${params.results_output}/results/mapped_bam/*mapped.realcells_only.bam")
        .map { file ->
            def sample = (file.name =~ /Isogut\d+/)[0]
            [sample, file]
        }
        .collect()
        .set { bam_files }
      bam_files.view()
}
