process SQANTI3_QC {
    label 'big_job'


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
