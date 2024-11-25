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
