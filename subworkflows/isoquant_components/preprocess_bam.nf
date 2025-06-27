
process preprocess_bam_perChr {

  label 'micro_multithread_job'

  input:
      tuple val(chrom), val(sample_id), path(bam), path(bai)
      val(suffix)
  output:
      tuple val(chrom), val(sample_id), path("${sample_id}.${chrom}${suffix}.bam"), path("${sample_id}.${chrom}${suffix}.bam.bai")
  script:
  """
  bash ${baseDir}/scripts/preprocess_bam.sh ${sample_id} ${bam} ${chrom} ${suffix} ${task.cpus}
  """

}

process preprocessWithExclusion_bam_perChr {

  label 'micro_multithread_job'

  input:
      tuple val(chrom), val(sample_id), path(bam), path(bai)
      val(suffix)
      path(bed)
  output:
      tuple val(chrom), val(sample_id), path("${sample_id}.${chrom}${suffix}.bam"), path("${sample_id}.${chrom}${suffix}.bam.bai")
  script:
  """
  samtools view -h -L ${bed} -U ${sample_id}.${chrom}${suffix}.tmp.bam ${bam} ${chrom} > /dev/null
  samtools index ${sample_id}.${chrom}${suffix}.tmp.bam

  bash ${baseDir}/scripts/preprocess_bam.sh ${sample_id} "${sample_id}.${chrom}${suffix}.tmp.bam" ${chrom} ${suffix} ${task.cpus}

  numreads_before_exclusion=\$(samtools view -c -h ${bam} ${chrom})
  numreads_after_exclusion=\$(samtools view -c -h ${sample_id}.${chrom}${suffix}.tmp.bam)
  numreads_after_processing=\$(samtools view -c -h ${sample_id}.${chrom}${suffix}.bam)

  echo "${chrom},${sample_id},\${numreads_before_exclusion},\${numreads_after_exclusion},\${numreads_after_processing}" > ${sample_id}_numreads.${chrom}.txt

  rm ${sample_id}.${chrom}${suffix}.tmp.bam
  rm ${sample_id}.${chrom}${suffix}.tmp.bam.bai
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
    if (params.isoquant_exclusion_regions_bed) {
      processedFullBamperChr_ch=preprocessWithExclusion_bam_perChr(preprocess_bam_perChr_input_ch,".mapped.realcells_only.processed",params.isoquant_exclusion_regions_bed)
    } else {
      processedFullBamperChr_ch=preprocess_bam_perChr(preprocess_bam_perChr_input_ch,".mapped.realcells_only.processed")
    }
  emit:
    processedFullBamperChr_ch

}
