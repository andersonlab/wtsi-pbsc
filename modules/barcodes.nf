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
