

workflow bamsWithExclusion {
  main:


    def excluded_samples_list = file(params.exclude_samples).exists() ?
      file(params.exclude_samples).text.split('\n').drop(1).collect { it.split(',')[0].trim() } : []




    ///Channel emitting sample_id,mapped BAM taking care to remove excluded samples
    Channel
      .fromPath(params.input_samples_path)
      .splitCsv(sep: ',', header: true)
      .filter { it -> !(it.sample_id in excluded_samples_list) }
      .map {it -> [it.sample_id,"${params.results_output}qc/mapped/"+it.sample_id+".mapped.realcells_only.bam","${params.results_output}qc/mapped/"+it.sample_id+".mapped.realcells_only.bam.bai"]}
      .set { fullBam_ch }
  emit:
    fullBam_ch

}
