process SQANTI3_QC {
    label 'big_job'

    publishDir "${params.results_output}results/transcript_info/sqanti3/", mode: 'copy', overwrite: true

    input:
      path(input_gtf_f)
      path(ref_gtf_f)
      path(genome_fasta_f)
      path(polya_f)
      path(cage_peak_f)
      path(polya_sites)
      val(sqanti3_path)

    output:
      path("sqanti3_qc")


    script:
      """
      python ${sqanti3_path}sqanti3_qc.py ${input_gtf_f} ${ref_gtf_f} ${genome_fasta_f} --polyA_motif_list ${polya_f} --CAGE_peak ${cage_peak_f} --report both -t ${task.cpus} --polyA_peak ${polya_sites} -d sqanti3_qc/ --force_id_ignore
      """

}

process SQANTI3_FILTER {
  label 'mini_job'

  publishDir "${params.results_output}results/transcript_info/sqanti3/", mode: 'copy', overwrite: true

  input:
    path(classification_f)
    path(input_gtf_f)
    val(sqanti3_path)
  output:
    path("sqanti3_filter")
  script:
  """
  python ${sqanti3_path}sqanti3_filter.py rules ${classification_f} -d sqanti3_filter/ -e
  python ${baseDir}/scripts/create_genedb.py -g ${input_gtf_f} -o sqanti3_filter/${input_gtf_f}.db
  python ${baseDir}/scripts/db_subset.py -d sqanti3_filter/${input_gtf_f}.db -i sqanti3_filter/${input_gtf_f.baseName}_inclusion-list.txt -o sqanti3_filter/${input_gtf_f.baseName}.filtered.gtf
  python ${baseDir}/scripts/create_genedb.py -g sqanti3_filter/${input_gtf_f.baseName}.filtered.gtf -o sqanti3_filter/${input_gtf_f.baseName}.filtered.gtf.db

  """

}
