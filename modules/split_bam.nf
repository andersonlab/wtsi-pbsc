process SPLIT_BAM_SINTO {
    label 'split_bam_sinto'

    input:
    tuple val(pool_id), path(bam), path(bai), path(barcode_list)

    output:
    tuple val(pool_id), path("split_bam/*.bam"), emit: donor_bams
    path "versions.yml", emit: versions

    script:
    """
    sinto filterbarcodes -b ${bam} -c ${barcode_list} -p ${task.cpus} --outdir split_bam
    all_barcode_groups=\$(cut -f 2 ${barcode_list} | sort | uniq)
    for barcode_group in \${all_barcode_groups}; do

        bash ${baseDir}/scripts/append_to_tag.sh -i split_bam/\${barcode_group}.bam -t CB:Z -s \${barcode_group} -O sam | bash ${baseDir}/scripts/append_to_readname.sh -s \${barcode_group} -O bam -o split_bam/\${barcode_group}.tmp.bam -n 

        mv split_bam/\${barcode_group}.tmp.bam split_bam/\${barcode_group}.bam
        samtools index -@ ${task.cpus} split_bam/\${barcode_group}.bam 

    done
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


