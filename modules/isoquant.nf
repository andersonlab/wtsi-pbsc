process isoquant_split_by_chr {
    label 'isoquant_split_by_chr'

    input:
        tuple val(sample_id), val(chrom), path(mapped_bam), path(mapped_bai)

    output:
    tuple val(sample_id), val(chrom), path("${sample_id}.${chrom}.mapped.realcells_only.bam"), path("${sample_id}.${chrom}.mapped.realcells_only.bam.bai"), emit: isoquant_split_tuple


    script:
    """
    samtools view -h ${mapped_bam} ${chrom} | awk -v sample=${sample_id} 'BEGIN {OFS="\t"}{if (\$1 ~ /^@/) {print; next}for(i=12; i<=NF; i++) {if (\$i ~ /^CB:Z:/) {\$i = \$i"_"sample}}print}' | samtools view -h -bo ${sample_id}.${chrom}.mapped.realcells_only.bam;
    samtools index ${sample_id}.${chrom}.mapped.realcells_only.bam
    """
}

process run_isoquant {
    label 'isoquant'

    input:
        tuple val(chrom), val(sample_ids), path(bams), path(bais)
        val gtf_f
        val genome_fasta_f

    output:
        tuple val(chrom), path("${chrom}/"), emit: isoquant_output

    script:
    """
    isoquant.py --reference ${genome_fasta_f} --genedb ${gtf_f} --complete_genedb --sqanti_output --bam ${bams.join(' ')} --labels ${sample_ids.join(' ')} --data_type pacbio_ccs -o ${chrom} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format linear --bam_tags CB --no_secondary --clean_start
    """
}
