
process preprocess_bam_perChr {

  label 'micro_multithread_job'

  input:
      tuple val(chrom), val(sample_id), path(bam), path(bai)
      val(suffix)
      val(bed_path)
  output:
      tuple val(chrom), val(sample_id), path("${sample_id}.${chrom}${suffix}.bam"), path("${sample_id}.${chrom}${suffix}.bam.bai")
  script:
  def has_bed = bed_path != ''
  def input_bam = has_bed ? "${sample_id}.${chrom}${suffix}.tmp.bam" : bam
  def pre_cmd = has_bed ? """
  samtools view -h -L ${bed_path} -U ${sample_id}.${chrom}${suffix}.tmp.bam ${bam} ${chrom} > /dev/null
  samtools index ${sample_id}.${chrom}${suffix}.tmp.bam
  """ : ""
  def post_cmd = has_bed ? """
  numreads_before_exclusion=\$(samtools view -c -h ${bam} ${chrom})
  numreads_after_exclusion=\$(samtools view -c -h ${sample_id}.${chrom}${suffix}.tmp.bam)
  numreads_after_processing=\$(samtools view -c -h ${sample_id}.${chrom}${suffix}.bam)
  echo "${chrom},${sample_id},\${numreads_before_exclusion},\${numreads_after_exclusion},\${numreads_after_processing}" > ${sample_id}_numreads.${chrom}.txt
  rm ${sample_id}.${chrom}${suffix}.tmp.bam
  rm ${sample_id}.${chrom}${suffix}.tmp.bam.bai
  """ : ""
  """
  ${pre_cmd}
  bash ${baseDir}/scripts/preprocess_bam.sh ${sample_id} ${input_bam} ${chrom} ${suffix} ${task.cpus}
  ${post_cmd}
  """

}


workflow preprocess_bam_perChr_wf {
  take:
    chrom_ch
    fullBam_ch
  main:
    chrom_ch
    .combine(fullBam_ch)
    .set {preprocess_bam_perChr_input_ch}
    def bed_path = params.isoquant_exclusion_regions_bed ?: ''
    processedFullBamperChr_ch = preprocess_bam_perChr(preprocess_bam_perChr_input_ch, ".mapped.realcells_only.processed", bed_path)
  emit:
    processedFullBamperChr_ch

}
