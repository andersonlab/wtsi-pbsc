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
