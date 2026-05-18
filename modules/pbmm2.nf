process PBMM2 {
    label 'map_reads'
    //publishDir "${params.results_output}qc/mapped", mode: 'copy'

    input:
        tuple val(sample_id), path(dedup_bam)
        val genome_fasta_f
        path manual_barcodes  // optional: pass [] if unused

    output:
        tuple val(sample_id), path("${dedup_bam.name.replaceAll(/\.bam/, '.mapped_chunk.realcells_only.bam')}"), emit: map_tuple

    script:
    def mapped_bam  = dedup_bam.name.replaceAll(/\.bam/, '.mapped_chunk.bam')
    def out_bam     = dedup_bam.name.replaceAll(/\.bam/, '.mapped_chunk.realcells_only.bam')
    def tag_manual  = (manual_barcodes) ? """
bash ${baseDir}/bin/tag_manual_barcodes.sh ${manual_barcodes} ${mapped_bam} ${mapped_bam}.tmp.bam
mv ${mapped_bam}.tmp.bam ${mapped_bam}
samtools index -@ ${task.cpus} ${mapped_bam}
""" : ""
    """
    samtools index -@ ${task.cpus} ${dedup_bam}
    pbmm2 align -j ${task.cpus} --preset ISOSEQ --sort ${dedup_bam} ${genome_fasta_f} ${mapped_bam}
    samtools index -@ ${task.cpus} ${mapped_bam}

    ${tag_manual}
    samtools view -@ ${task.cpus} -h -d rc:1 -bo ${out_bam} ${mapped_bam}
    """
}
