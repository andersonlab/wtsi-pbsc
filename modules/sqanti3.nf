process SQANTI3_QC_SPLIT {
    label 'sqanti3_qc_split'

    input:
      path(input_gtf_f)
      val(num_chunks)

    output:
      path("chunks/chunk_*.gtf"), emit: chunk_gtfs

    script:
      """
        python ${baseDir}/scripts/sqanti3_gtf_chunks.py -g ${input_gtf_f} -n ${num_chunks} -o chunks/
      """
}

process SQANTI3_QC_CHUNK {
    label 'sqanti3_qc_chunk'

    input:
      path(chunk_gtf)
      path(ref_gtf_f)
      path(genome_fasta_f)
      path(polya_f)
      path(cage_peak_f)
      path(polya_sites)
      val(num_chunks)

    output:
      path("qc_chunk_*"), emit: chunk_dir

    script:
      """
        CHUNK_IDX=\$(basename ${chunk_gtf} .gtf | sed 's/^chunk_//')
        mkdir -p qc_chunk_\${CHUNK_IDX}
        export OMP_NUM_THREADS=${task.cpus}
        export MKL_NUM_THREADS=${task.cpus}
        python ${baseDir}/scripts/sqanti3_qc_chunk_run.py \
              --isoforms ${chunk_gtf} --refGTF ${ref_gtf_f} --refFasta ${genome_fasta_f} \
              --polyA_motif_list ${polya_f} --CAGE_peak ${cage_peak_f} --polyA_peak ${polya_sites} \
              --chunks ${num_chunks} -t ${task.cpus} \
              -d qc_chunk_\${CHUNK_IDX}/ --include_ORF --output transcript_models --report skip
      """
}

process SQANTI3_QC_COMBINE {
    label 'sqanti3_qc_combine'

    publishDir "${params.results_output}results/transcript_info/sqanti3/", mode: 'copy', overwrite: true

    input:
      path(chunk_dirs)
      path(input_gtf_f)
      path(ref_gtf_f)
      path(genome_fasta_f)
      val(num_chunks)

    output:
      path("sqanti3_qc"),                                        emit: qc_dir
      path("sqanti3_qc/transcript_models_classification.txt"),   emit: classification
      path("sqanti3_qc/transcript_models_corrected.gtf"),        emit: corrected_gtf
      path("sqanti3_qc/transcript_models_corrected.fasta"),      emit: corrected_fasta
      path("sqanti3_qc/transcript_models_corrected.faa"),        emit: corrected_faa
      path("sqanti3_qc/TD2"),                                    emit: td2_dir

    script:
      """
        python ${baseDir}/scripts/sqanti3_qc_combine.py \
              --isoforms ${input_gtf_f} --refGTF ${ref_gtf_f} --refFasta ${genome_fasta_f} \
              --chunks ${num_chunks} -t ${task.cpus} \
              -d sqanti3_qc/ --include_ORF --output transcript_models --report pdf \
              --chunk_dirs ${chunk_dirs.join(' ')}
      """
}

process SQANTI3_FILTER {
  label 'sqanti3_filter'

  publishDir "${params.results_output}results/transcript_info/sqanti3/", mode: 'copy', overwrite: true

  input:
    path(classification_f)
    path(corrected_gtf_f)
    path(corrected_faa)
    path(td2_dir)
    path(input_gtf_f)
    path(sqanti_filter_json)
  output:
    path("sqanti3_filter"),                                                  emit: filter_dir
    path("sqanti3_filter/transcript_models_pass_isoforms.txt"),              emit: pass_isoforms
    path("sqanti3_filter/transcript_models_corrected.filtered.gtf"),         emit: filtered_gtf
    path("sqanti3_filter/transcript_models.filtered.gtf"),                   emit: filtered_original_gtf
    path("sqanti3_filter/TD2"),                                              emit: filtered_td2_dir
  script:
  def prefix = classification_f.baseName.replace("_classification", "")
  """
  sqanti3_filter.py rules --sqanti_class ${classification_f} \
      -j ${sqanti_filter_json} -d sqanti3_filter/ --filter_gtf ${corrected_gtf_f} --filter_faa ${corrected_faa} --cpus ${task.cpus} --skip_report
  python ${baseDir}/scripts/gtf_subset.py -g ${corrected_gtf_f} -i sqanti3_filter/transcript_models_pass_isoforms.txt -o sqanti3_filter/${prefix}_corrected.filtered.gtf
  python ${baseDir}/scripts/gtf_subset.py -g ${input_gtf_f} -i sqanti3_filter/transcript_models_pass_isoforms.txt -o sqanti3_filter/${prefix}.filtered.gtf
  python ${baseDir}/scripts/td2_filter.py -d ${td2_dir} -i sqanti3_filter/transcript_models_pass_isoforms.txt -o sqanti3_filter/TD2
  """

}
