process MTX_SUBSET {
    label 'small_job'


    input:
      path(input_mtx_dir)
      path(subset_f)

    output:
      path("${input_mtx_dir}_filtered")
    script:
      """
      mkdir -p ${input_mtx_dir}_filtered
      python ${baseDir}/scripts/mtx_subset.py -i ${input_mtx_dir}/ -s ${subset_f} -d ${input_mtx_dir}_filtered/
      """

}
