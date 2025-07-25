process barcode_correction {
    label 'barcodes'
    publishDir "${params.results_output}qc/correct", mode: 'copy'

    input:
        tuple val(sample_id), path(refined_reads_bam)
        path(threeprime_whitelist)
        val barcode_correction_method
        val barcode_correction_percentile

    output:
        tuple val(sample_id), path("${sample_id}.corrected.sorted.bam"), emit: barcode_corrected_tuple
        path "*corrected.sorted*"
        path "*corrected.report*"


    script:
        """
        # Correct step
        if [[ "${barcode_correction_method}" == "percentile" ]]; then
            isoseq correct  -j ${task.cpus} --method "${barcode_correction_method}" --percentile "${barcode_correction_percentile}" --barcodes "${threeprime_whitelist}" "${refined_reads_bam}" "${sample_id}.corrected.bam"
        elif [[ "${barcode_correction_method}" == "knee" ]]; then
            isoseq correct -j ${task.cpus} --method "${barcode_correction_method}" --barcodes "${threeprime_whitelist}" "${refined_reads_bam}" "${sample_id}.corrected.bam"
        else
            echo "Invalid barcode correction method: ${barcode_correction_method}" >&2
            exit 1
        fi

        # Sort step
        samtools sort -@ ${task.cpus} -t CB "${sample_id}.corrected.bam" -o "${sample_id}.corrected.sorted.bam"
        samtools index "${sample_id}.corrected.sorted.bam"

        # BCStats step
        if [[ "${barcode_correction_method}" == "percentile" ]]; then
            isoseq bcstats --method "${barcode_correction_method}" --percentile "${barcode_correction_percentile}" --json "${sample_id}.corrected.sorted.json" -o "${sample_id}.corrected.sorted.tsv" "${sample_id}.corrected.sorted.bam"
        elif [[ "${barcode_correction_method}" == "knee" ]]; then
            isoseq bcstats --method "${barcode_correction_method}" --json "${sample_id}.corrected.sorted.json" -o "${sample_id}.corrected.sorted.tsv" "${sample_id}.corrected.sorted.bam"
        else
            echo "Invalid barcode correction method: ${barcode_correction_method}" >&2
            exit 1
        fi
        """

}

process get_barcodes {
    label 'process_low'
    
    input:
    tuple val(sample), path(bam)
    val(N)

    output:
    tuple val(sample), path("barcodes_*.txt"), emit: barcodes_tuple

    script:
        def K = N.toString().length() + 1
        """ 
        	samtools view ${bam} | awk '{ for(i=12;i<=NF;i++) if(\$i ~ /^CB:Z:/) {split(\$i,a,":"); print a[3]} }' | sort | uniq > tmp.txt
            #samtools view ${bam} | grep -o 'CB:Z:[^[:space:]]*' | cut -d: -f3 | sort | uniq > tmp.txt
            split -n l/${N} -d -a ${K} --additional-suffix=.txt tmp.txt barcodes_
        """

}


process supset_bam {
    label 'supset_bam'
       
    input:
        tuple val(sample), path(bam), path(barcodes)

    output:
        tuple val(sample), path("*.splited.bam"), emit: chunk_tuple

    script:
    	def bam_name=bam.baseName
    	def barcode_name=barcodes.baseName
        """ 
            samtools view --threads ${task.cpus} --tag-file CB:${barcodes} \
               -o ${bam_name}.${barcode_name}.splited.bam ${bam}
            samtools index -c ${bam_name}.${barcode_name}.splited.bam
        """

}

process supset_bam_with_bai {
    label 'supset_bam'
       
    input:
        tuple val(sample), path(bam), path(bai), path(barcodes)

    output:
        tuple val({ "${sample}___${barcodes.simpleName.replaceFirst('barcodes__','')}" }), path("*.splited.bam"), path("*.splited.bam.bai"), emit: per_donor_tuple

    script:
    	def bam_name=bam.baseName
    	def barcode_name=barcodes.baseName.replaceAll('barcodes__','')
        """ 
            samtools view --threads ${task.cpus} --tag-file CB:${barcodes} \
               -o ${bam_name}.${barcode_name}.splited.bam ${bam}
            samtools index ${bam_name}.${barcode_name}.splited.bam
        """

}

process dedup_reads {
    label 'deduplication'

    input:
        tuple val(sample_id), path(barcode_corrected_chunk_bam)
        val (b_size)


    output:
        tuple val(sample_id), path("${barcode_corrected_chunk_bam.name.replaceAll(/\.bam/, '.dedup.bam')}"), emit: dedup_tuple

    script:
    """
    isoseq groupdedup  -j ${task.cpus} --batch-size ${b_size} --keep-non-real-cells ${barcode_corrected_chunk_bam} ${barcode_corrected_chunk_bam.name.replaceAll(/\.bam/, '.dedup.bam')}

    """
}


process combine_dedups {
    label 'combine_bams'
    publishDir "${params.results_output}qc/dedup", mode: 'copy'

    input:
        tuple val(sample_id), path(dedup_bam_chunks)


    output:
        tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: dedup_tuple

    script:
    """
 	samtools merge -f ${sample_id}.dedup.bam ${dedup_bam_chunks.join(' ')}
    samtools index -@ ${task.cpus} ${sample_id}.dedup.bam

    """
}

process combine_mupped {
    label 'combine_bams'
    publishDir "${params.results_output}qc/mapped", mode: 'copy'

    input:
        tuple val(sample_id), path(mapped_bam_chunks)


    output:
        tuple val(sample_id), path("${sample_id}.mapped.realcells_only.bam"), path("${sample_id}.mapped.realcells_only.bam.bai"), emit: dedup_tuple

    script:
    """
 	samtools merge -f ${sample_id}.mapped.realcells_only.merged.bam ${mapped_bam_chunks.join(' ')}
    samtools sort -@ ${task.cpus} -o ${sample_id}.mapped.realcells_only.bam ${sample_id}.mapped.realcells_only.merged.bam
    samtools index -@ ${task.cpus} ${sample_id}.mapped.realcells_only.bam
    """
}


process bam_stats {
    label 'bam_stats'
    publishDir "${params.results_output}qc/dedup", mode: 'copy'

    input:
        tuple val(sample_id), path(dedup_bam), path("${sample_id}.dedup.bam.bai")
        val barcode_correction_method
        val barcode_correction_percentile


    output:
        tuple val(sample_id), path(dedup_bam), path("${sample_id}.dedup.bam.bai")

    script:
    """
    if [[ "${barcode_correction_method}" == "percentile" ]]; then
      isoseq bcstats  -j ${task.cpus} --method ${barcode_correction_method} --percentile ${barcode_correction_percentile} --json ${sample_id}.dedup.json -o ${sample_id}.dedup.tsv ${sample_id}.dedup.bam
    elif [[ "${barcode_correction_method}" == "knee" ]]; then
      isoseq bcstats  -j ${task.cpus} --method ${barcode_correction_method} --json ${sample_id}.dedup.json -o ${sample_id}.dedup.tsv ${sample_id}.dedup.bam
    else
        echo "Invalid barcode correction method: ${barcode_correction_method}" >&2
        exit 1
    fi;

    """
}


