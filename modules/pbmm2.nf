process pbmm2 {
    label 'map_reads'
    publishDir "${params.results_output}qc/mapped", mode: 'copy'

    input:
        tuple val(sample_id), path(dedup_bam), path(bam_bai)
        val genome_fasta_f

    output:
        tuple val(sample_id), path("${sample_id}.mapped.realcells_only.bam"), path("${sample_id}.mapped.realcells_only.bam.bai")

    script:
    """
    pbmm2 align -j ${task.cpus} --preset ISOSEQ --sort ${dedup_bam} ${genome_fasta_f} ${sample_id}.mapped.bam
    samtools index -@ ${task.cpus} ${sample_id}.mapped.bam

    samtools view -@ ${task.cpus} -h -d rc:1 -bo ${sample_id}.mapped.realcells_only.bam ${sample_id}.mapped.bam
    samtools index -@ ${task.cpus} ${sample_id}.mapped.realcells_only.bam
    """
}
