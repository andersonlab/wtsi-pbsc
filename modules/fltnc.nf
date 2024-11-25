process split_reads {

    label 'split_reads'

    publishDir "${params.results_output}qc/skera_reports", mode: 'copy', pattern: '*.segmented.summary.json'

    input:
	   tuple val(sample_id), path(input_bam)
     val skera_primers

    output:
      tuple val(sample_id), path("${sample_id}.segmented.bam"), emit: split_reads_tuple
    script:
    """
    pbskera split ${input_bam} ${skera_primers} ${sample_id}.segmented.bam
    """

    stub:
    """
    touch ${sample_id}.segmented.bam
    pbskera split --help
    """
}

process remove_primer {
    label 'remove_primer'

    publishDir "${params.results_output}qc/lima_reports", mode: 'copy', pattern: "*.lima.summary"

    input:
      tuple val(sample_id), path(segmented_bam)
      val tenx_primers

    output:
      tuple val(sample_id), path("${sample_id}.5p--3p.bam"), emit: removed_primer_tuple

    script:
    """
    lima ${segmented_bam} ${tenx_primers} ${sample_id}.bam --isoseq
    """
    stub:
    """
    touch ${sample_id}.5p--3p.bam
    lima --help
    """


}

process tag_bam {
    label 'remove_primer'

    input:
      tuple val(sample_id), path(primer_removed_bam)

    output:
    tuple val(sample_id), path("${sample_id}.flt.bam"), emit: tagged_tuple

    script:
    """
    isoseq tag ${primer_removed_bam} ${sample_id}.flt.bam --design T-12U-16B
    """
    stub:
    """
    touch ${sample_id}.flt.bam
    isoseq tag --help
    """

}

process refine_reads {
    label 'remove_primer'

    publishDir "${params.results_output}qc/refined", mode: 'copy'


    input:
      tuple val(sample_id), path(tagged_bam)
      val tenx_primers
      val min_polya_length

    output:
      path "${sample_id}.fltnc.bam"
      path "${sample_id}.fltnc.bai"
      path "${sample_id}.fltnc.filter_summary.report.json"

    script:
    """
    isoseq refine ${tagged_bam} ${tenx_primers} ${sample_id}.fltnc.bam --require-polya --min-polya-length ${min_polya_length}
    samtools index ${sample_id}.fltnc.bam
    """
    stub:
    """
    touch ${sample_id}.fltnc.bam
    """
}
