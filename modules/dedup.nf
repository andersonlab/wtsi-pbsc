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
    echo "test"
    isoseq groupdedup -j ${task.cpus} --keep-non-real-cells ${sorted_barcodes_bam} \$sample_name.dedup.bam
    isoseq bcstats -j ${task.cpus} --method ${params.barcode_correction_method} --json \$sample_name.dedup.json -o \$sample_name.dedup.tsv ${dedup_bam}
    samtools index \$sample_name.dedup.bam
    """
}
