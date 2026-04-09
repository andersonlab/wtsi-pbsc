//process GET_BARCODES{
//  label 'get_barcodes'
//
//  input:
//    tuple val(pool_id), path(vireo_tsv)
//
//  output:
//    tuple val(pool_id), path("*bc_list.txt.gz"), emit: barcode_lists
//    //path "versions.yml", emit: versions
//
//  script:
//  """
//    awk -F'\t' -v pool="${pool_id}" \
//    'NR>1 && \$2 ~ /^donor[0-9]+\$/ && !seen[\$2,\$1]++ {
//        file = pool "_" \$2 ".bc_list.txt"
//        print \$1 > file
//    }' ${vireo_tsv}
//    gzip ${pool_id}_donor*.bc_list.txt
//
//
//    #cat <<-END_VERSIONS > versions.yml
//    #"${task.process}":
//    #    python: \$(python --version | sed 's/Python //g')
//    #END_VERSIONS
//  """
//}


process GET_BARCODES{
  label 'get_barcodes'

  input:
    tuple val(pool_id), path(vireo_tsv)

  output:
    tuple val(pool_id), path("${pool_id}.bc_list.txt.gz"), emit: barcode_list
    //path "versions.yml", emit: versions

  script:
  """

    tail -n +2 ${vireo_tsv} | awk -F'\t' -v pool="${pool_id}" 'BEGIN {OFS="\t"} \$2 !~ /^(unassigned|doublet)\$/ {gsub(/^donor/, pool"_donor", \$2); print \$1, \$2}' > ${pool_id}.bc_list.txt
    gzip ${pool_id}.bc_list.txt

    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    python: \$(python --version | sed 's/Python //g')
    #END_VERSIONS
  """
}



process SPLIT_BAM {
    label 'split_bam'

    input:
    tuple val(pool_id), path(bam), path(bai), path(barcode_list)

    output:
    tuple val(pool_id), path("split_bam/*.bam"), emit: donor_bams
    path "versions.yml", emit: versions

    script:
    """
    sinto filterbarcodes -b ${bam} -c ${barcode_list} -p ${task.cpus} --outdir split_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sinto: \$(sinto --version)
    END_VERSIONS
    """
}


process INDEX_SPLIT_BAM {
    label 'index_bam'
    publishDir  path: "${params.results_output}deconvolution/bam/",
          pattern: "*.bam",
          mode: 'copy',
          overwrite: "true"
    publishDir  path: "${params.results_output}deconvolution/bam/",
          pattern: "*.bai",
          mode: 'copy',
          overwrite: "true"

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path(bam), path("${bam}.bai"), emit: donor_bams
    path "versions.yml", emit: versions

    script:
    """
    samtools index -@ ${task.cpus} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1)
    END_VERSIONS
    """
}


