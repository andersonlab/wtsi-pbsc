process pbmm2 {
    label 'map_reads'
    //publishDir "${params.results_output}qc/mapped", mode: 'copy'

    input:
        tuple val(sample_id), path(dedup_bam)
        val genome_fasta_f

    output:
        tuple val(sample_id), path("${dedup_bam.name.replaceAll(/\.bam/, '.mapped_chunk.realcells_only.bam')}"), emit: map_tuple

    script:
    """
    samtools index -@ ${task.cpus} ${dedup_bam}
    pbmm2 align -j ${task.cpus} --preset ISOSEQ --sort ${dedup_bam} ${genome_fasta_f} ${dedup_bam.name.replaceAll(/\.bam/, '.mapped_chunk.bam')}
    samtools index -@ ${task.cpus} ${dedup_bam.name.replaceAll(/\.bam/, '.mapped_chunk.bam')}

    samtools view -@ ${task.cpus} -h -d rc:1 -bo ${dedup_bam.name.replaceAll(/\.bam/, '.mapped_chunk.realcells_only.bam')} ${dedup_bam.name.replaceAll(/\.bam/, '.mapped_chunk.bam')}
    """
}
