nextflow.enable.dsl=2

/////// Parameters ///////////
params.input_samples_path = '/nfs/team152/oe2/isogut/scripts/workflows/isoseq/input_samples.tsv'
params.skera_primers='/lustre/scratch126/humgen/projects/isogut/qc/mas16_primers.fasta'
params.tenx_primers='/lustre/scratch126/humgen/projects/isogut/qc/3_prime_kit_primers.fasta'
params.results_output="/lustre/scratch126/humgen/projects/isogut/output/"
params.threeprime_whitelist="/lustre/scratch126/humgen/projects/isogut/qc/10x_whitelist_3_prime"
params.genome_fasta_f="/lustre/scratch126/humgen/projects/isogut/qc/human_GRCh38_no_alt_analysis_set.fasta"
params.cage_peak_f="/lustre/scratch126/humgen/projects/isogut/qc/refTSS_v3.3_human_coordinate.hg38.sorted.bed"
params.polya_f="/lustre/scratch126/humgen/projects/isogut/qc/polyA.list.txt"
params.barcode_correction_method="percentile"
params.gtf_f="/lustre/scratch126/humgen/projects/isogut/utils/gencode.v46.annotation.sorted.gtf"
params.sqanti3_path="/nfs/team152/oe2/isogut/software/SQANTI3-5.2.1/"
params.polya_sites="/lustre/scratch126/humgen/projects/isogut/utils/atlas.clusters.2.0.GRCh38.96.bed"
params.utils_scripts_path="/nfs/team152/oe2/isogut/scripts/utils/"


params.knee_plot_script='/nfs/team152/oe2/isogut/scripts/utils/plot_knees.py'


#This will replace the current hard-coded params
if (params.config) {
    includeConfig(params.config)
}
//////////////////////////////
process split_reads {

    label 'skera'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}qc/split_reads", mode: 'copy'

    input:
        path bam_path

    output:
        path "*.segmented.bam", emit: split_reads_bam
        path "*.segmented.bam.bai"
        path "*.segmented.summary.csv"
        path "*.segmented.read_lengths.csv"
    script:
    """
    sample_name=\$(echo ${bam_path} | cut -d'.' -f1)
    skera split ${bam_path} ${params.skera_primers}  \$sample_name.segmented.bam --num-threads 20
    samtools index \$sample_name.segmented.bam
    """


}

process remove_primer {
    label 'remove_primer'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}qc/removed_primer", mode: 'copy'

    input:
        path segmented_bam_path

    output:
        path "*.5p--3p.bam", emit: removed_primer_bam
        path "*.5p--3p.bam.bai"
        path "*.lima.report"
        path "*.lima.summary"

    script:
    """
    sample_name=\$(echo ${segmented_bam_path} | cut -d'.' -f1)
    lima ${segmented_bam_path} ${params.tenx_primers} \$sample_name.bam --isoseq --num-threads 20
    samtools index \$sample_name.5p--3p.bam
    """


}

process remove_tags {
    label 'remove_primer'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}qc/tagged", mode: 'copy'

    input:
        path primer_removed_bam_path

    output:
    path "*.flt.bam", emit: removed_tag_bam
    path "*.flt.bam.bai"

    script:
    """
    sample_name=\$(echo ${primer_removed_bam_path} | cut -d'.' -f1)
    isoseq tag ${primer_removed_bam_path} \$sample_name.flt.bam --design T-12U-16B
    samtools index \$sample_name.flt.bam
    """

}

process refine_reads {
    label 'remove_primer'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}qc/refined", mode: 'copy'


    input:
        path removed_tag_bam

    output:
        path "*fltnc.bam", emit: refined_bam
        path "*fltnc.bam.bai"
        path "*fltnc.bam.pbi"
        path "*fltnc.report.csv"
        path "*fltnc.filter_summary.report.json"
        path "*fltnc.consensusreadset.xml"

    script:
    """
    sample_name=\$(echo ${removed_tag_bam} | cut -d'.' -f1)
    isoseq refine ${removed_tag_bam} ${params.tenx_primers} \$sample_name.fltnc.bam --require-polya
    samtools index \$sample_name.fltnc.bam
    rm ${removed_tag_bam}
    """
}

process refine_reads__withpolya {
    label 'remove_primer'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}qc/refined_withpolya", mode: 'copy'


    input:
        path removed_tag_bam

    output:
        path "*fltnc.bam", emit: refined_bam
        path "*"

    script:
    """
    sample_name=\$(echo ${removed_tag_bam} | cut -d'.' -f1)
    isoseq refine ${removed_tag_bam} ${params.tenx_primers} \$sample_name.fltnc.bam --require-polya
    samtools index \$sample_name.fltnc.bam
    rm ${removed_tag_bam}
    """
}

process postrefine_stats {
    label 'bcstats'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}qc/refined", mode: 'copy'

    input:
        path refined_bam

    output:
        path "*fltnc.knee.json"
        path "*fltnc.knee.tsv"
        path "*fltnc.knee.knee.png"
        path "*fltnc.perc.json"
        path "*fltnc.perc.tsv"
        path "*fltnc.perc.knee.png"

    script:
    """
    sample_name=\$(echo ${refined_bam} | cut -d'.' -f1)
    isoseq bcstats --json \${sample_name}.fltnc.perc.json -o \${sample_name}.fltnc.perc.tsv --method percentile ${refined_bam}
    isoseq bcstats --json \${sample_name}.fltnc.knee.json -o \${sample_name}.fltnc.knee.tsv --method knee ${refined_bam}

    python ${params.knee_plot_script} -t \${sample_name}.fltnc.perc.tsv -o \${sample_name}.fltnc.perc --estimate_percentile 95
    python ${params.knee_plot_script} -t \${sample_name}.fltnc.knee.tsv -o \${sample_name}.fltnc.knee
    """
}

process barcode_correction {
    label 'barcodes'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    input:
        path refined_reads_bam

    output:
        path "*corrected.bam", emit: barcode_corrected_bam
        path "*"


    script:
    """
    sample_name=\$(echo ${refined_reads_bam} | cut -d'.' -f1)
    isoseq correct --method ${params.barcode_correction_method} --barcodes ${params.threeprime_whitelist} ${refined_reads_bam} \$sample_name.corrected.bam
    rm ${refined_reads_bam}
    """
}

process sort_cellbarcodes {
    label 'sort_bc'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
    publishDir "${params.results_output}qc/corrected", mode: 'copy'

    input:
        path barcode_corrected_bam

    output:
        path "*corrected.sorted.bam", emit: sorted_bam
        path "*corrected.sorted.json"
        path "*corrected.sorted.tsv"
        path "*corrected.sorted.bam.bai"


    script:
    """
    sample_name=\$(echo ${barcode_corrected_bam} | cut -d'.' -f1)
    samtools sort -t CB ${barcode_corrected_bam} -o \$sample_name.corrected.sorted.bam
    isoseq bcstats --method ${params.barcode_correction_method} --json \$sample_name.corrected.sorted.json -o \$sample_name.corrected.sorted.tsv \$sample_name.corrected.sorted.bam
    samtools index \$sample_name.corrected.sorted.bam
    rm ${barcode_corrected_bam}
    """
}

process dedup_reads {
    label 'deduplication'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}results/unmapped_bam", mode: 'copy'

    input:
        path sorted_barcodes_bam

    output:
        path "*dedup.bam", emit: dedup_bam
        path "*dedup.bam.bai"
        path "*dedup.bam.pbi"
        path "*dedup.fasta"

    script:
    """
    sample_name=\$(echo ${sorted_barcodes_bam} | cut -d'.' -f1)
    isoseq groupdedup --keep-non-real-cells ${sorted_barcodes_bam} \$sample_name.dedup.bam
    samtools index \$sample_name.dedup.bam
    rm ${sorted_barcodes_bam}
    """
}

process postdedup_stats {
    label 'remove_primer'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
    publishDir "${params.results_output}results/unmapped_bam", mode: 'copy'

    input:
        path dedup_bam

    output:
        path "*dedup.json"
        path "*dedup.tsv"

    script:
    """
    sample_name=\$(echo ${dedup_bam} | cut -d'.' -f1)
    isoseq bcstats --method ${params.barcode_correction_method} --json \$sample_name.dedup.json -o \$sample_name.dedup.tsv ${dedup_bam}
    """
}




// Processes for the isoform classification workflow
process map_reads {
    label 'map_reads'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}results/mapped_bam", mode: 'copy'

    input:
        path unmapped_bam

    output:
        path "*mapped.bam", emit: mapped_bam
        path "*mapped.bam.bai"

    script:
    """
    sample_name=\$(echo ${unmapped_bam} | cut -d'.' -f1)
    pbmm2 align --preset ISOSEQ --sort ${unmapped_bam} ${params.genome_fasta_f} \$sample_name.mapped.bam
    samtools index \$sample_name.mapped.bam
    rm ${unmapped_bam}
    """
}
process make_real_cells_bams {
  label 'small_job'

  conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
  publishDir "${params.results_output}results/mapped_bam", mode: 'copy'

  input:
      path mapped_bam

  output:
      path "*mapped.realcells_only.bam", emit: mapped_realcells_only_bam
      path "*mapped.realcells_only.bam.bai"
  script:
  """
  sample_name=\$(echo ${mapped_bam} | cut -d'.' -f1)
  samtools view -h -d rc:1 -bo \$sample_name.mapped.realcells_only.bam ${mapped_bam}
  samtools index \$sample_name.mapped.realcells_only.bam
  """

}
process create_bam_file_list {
    label 'small_job'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
    publishDir "${params.results_output}results/mapped_bam", mode: 'copy'

    input:
    val realcells_bam_list

    output:
    path 'realcells_bam_list.fofn'

    script:
    """
    echo '${realcells_bam_list.join("\n")}' > realcells_bam_list.fofn
    """
}

process create_refined_bam_file_list {
    label 'small_job'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
    publishDir "${params.results_output}qc/refined", mode: 'copy'

    input:
    val bam_list

    output:
    path 'refined_bam_list.fofn'

    script:
    """
    echo '${bam_list.join("\n")}' > refined_bam_list.fofn
    """
}
process cluster_bams {
  label 'cluster_bams'
  conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
  publishDir "${params.results_output}results/clustered", mode: 'copy'


  input:
  path bam_file_list

  output:
  path "clustered.bam", emit: clustered_bam
  path "clustered.bam.pbi"

  script:
  """
  isoseq cluster2 ${bam_file_list} clustered.bam
  """

}
process collapse_isoforms {
    label 'map_reads'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
    publishDir "${params.results_output}results/collapsed", mode: 'copy'

    input:
        path mapped_reads_path

    output:
        path "*collapsed.gff", emit: 'collapsed_gff_path'
        path "*collapsed.group.txt", emit: 'collapsed_group_path'
        path "*collapsed.abundance.txt", emit: 'collapsed_abundance_path'

    script:
    """
    sample_name=\$(echo ${mapped_reads_path} | cut -d'.' -f1)
    isoseq collapse ${mapped_reads_path} \$sample_name.collapsed.gff
    rm ${mapped_reads_path}
    """
}

process prepare_input_GFF {
    label 'skera'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
    publishDir "${params.results_output}results/prepared", mode: 'copy'

    input:
        path collapsed_isoforms_path

    output:
        path "*collapsed.sorted.gff"

    script:
    """
    pigeon prepare ${collapsed_isoforms_path}
    rm ${collapsed_isoforms_path}
    """
}
process sqanti3_classify_isoforms {


  label 'sqanti3'

  conda "/software/hgi/envs/conda/team152/oe2/SQANTI3.env"

  publishDir { "${params.results_output}results/sqanti3_qc/${sample_name}" }, mode: 'copy'


  input:
      tuple val(sample_name), path(collapsed_sorted_isoform_path), path(collapsed_abundance_path)

    output:
        path '*_classification.txt', emit: 'classification_files'
        path '*_junctions.txt', emit: 'junction_files'
        path '*_corrected.faa', emit: 'faa_files'
        path '*_corrected.gtf', emit: gtf_files
        path '*'


  script:
  """
  python ${params.sqanti3_path}sqanti3_qc.py ${collapsed_sorted_isoform_path} ${params.gtf_f} ${params.genome_fasta_f} -fl ${collapsed_abundance_path} -o ${sample_name} --polyA_motif_list ${params.polya_f} --CAGE_peak ${params.cage_peak_f} --report pdf -t ${task.cpus} --polyA_peak ${params.polya_sites}
  """
}

process sqanti3_filter_isoforms {


  label 'sqanti3'

  conda "/software/hgi/envs/conda/team152/oe2/SQANTI3.env"

  publishDir { "${params.results_output}results/sqanti3_filter/${sample_name}" }, mode: 'copy'


  input:
      tuple val(sample_name), path(classification_f), path(sample_gtf_f), path(sample_faa_f)

    output:
        path '*_RulesFilter_result_classification.txt', emit: 'filtered_classification_files'
        path '*_inclusion-list.txt', emit: 'filtered_inclusion_list_files'
        path '*'


  script:
  """
  python ${params.sqanti3_path}sqanti3_filter.py rules ${classification_f} -o ${sample_name} --gtf ${sample_gtf_f} --faa ${sample_faa_f} -e -j ${params.utils_scripts_path}sqanti3_filter.json
  """
}

process sqanti3_seurat_isoform {
    label 'big_job'

    conda "/software/hgi/envs/conda/team152/oe2/isoseq/"

    publishDir "${params.results_output}results/sqanti3_seurat", mode: 'copy'

    input:
        tuple val(sample_name), path(fasta), path(unfiltred_classification), path(filtered_classification), path(abundance), path(group)

    output:
        path "*_RulesFilter_result_classification_pigeonformatted.txt", emit: pigeonformatted_filtered_classification_files
        path "*"

    script:
    """
    python ${params.utils_scripts_path}merge_abundance_sqanti3.py -c ${unfiltred_classification} -a ${abundance} -o ${sample_name}
    pigeonformatted_unfiltered_classification_f=${sample_name}_unfiltered_classification_pigeonformatted.txt
    python ${params.utils_scripts_path}merge_pigeonformatted_classification_filter_files.py -c \${pigeonformatted_unfiltered_classification_f} -f ${filtered_classification} -o ${sample_name}
    pigeonformatted_unfiltered_classification_f=${sample_name}_RulesFilter_result_classification_pigeonformatted.txt
    pigeon make-seurat --dedup ${fasta} --group ${group} \${pigeonformatted_unfiltered_classification_f} --keep-ribo-mito-genes -d ${sample_name} -o ${sample_name}
    """
}

process classify_isoforms {
    label 'classify_isoforms'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
    publishDir "${params.results_output}results/classified", mode: 'copy'

    input:
        path collapsed_sorted_isoform_path
        path collapsed_abundance_path

    output:
        path "*_classification.txt", emit: 'classification_files'
        path "*_junctions.txt", emit: 'junction_files'
        path "*summary.txt"
        path "*report.json"

    script:
    """
    pigeon classify ${collapsed_sorted_isoform_path} ${params.gtf_f} ${params.genome_fasta_f} --cage-peak ${params.cage_peak_f} --poly-a ${params.polya_f} --flnc ${collapsed_abundance_path} --gene-id --log-level DEBUG
    rm ${collapsed_sorted_isoform_path}
    """
}

process filter_isoforms {
    label 'classify_isoforms'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"
    publishDir "${params.results_output}results/filtered", mode: 'copy'

    input:
        path classification_files
        path junction_files
        path collapsed_sorted_isoforms_path

    output:
        path "*.filtered_lite_classification.txt", emit: 'classification_filtered_files'
        path "*.filtered_lite_junctions.txt"
        path "*.filtered_lite_reasons.txt"
        path "*.filtered.report.json"
        path "*.filtered.summary.txt"

    script:
    """

    pigeon filter ${classification_files} ${junction_files} --isoforms ${collapsed_sorted_isoforms_path}
    """
}

process seurat_isoform {
    label 'seurat_isoforms'

    conda "/software/hgi/envs/conda/team152/mt27/isoseq/"

    publishDir "${params.results_output}results/seurat", mode: 'copy'

    input:
        tuple val(sample_name), path(fasta), path(classification), path(group)

    output:
        path "*"

    script:
    """
    sample_name=\$(echo ${fasta} | cut -d'.' -f1)
    pigeon make-seurat --dedup ${fasta} --group ${group} ${classification} --keep-ribo-mito-genes -d \$sample_name
    """
}






/// This QC workflow will do the pre-processing in: https://isoseq.how/umi/cli-workflow.html

workflow split_tag_refine {

    Channel
        .fromPath(params.input_samples_path)
        .splitCsv(sep: ',' , header: true)
        .map { it -> it.long_read_path} // Extract the long_read_path column
        .set {bam_path}

    reads_split       = split_reads(bam_path)
    primer_removed    = remove_primer(reads_split.split_reads_bam)
    tag_removed       = remove_tags(primer_removed.removed_primer_bam)
    refined_reads     = refine_reads(tag_removed.removed_tag_bam)
    refined_bam_stats = postrefine_stats(refined_reads.refined_bam)
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
