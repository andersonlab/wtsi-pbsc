process barcode_correction {
    label 'barcodes'

    input:
        tuple val(sample_id), path(refined_reads_bam)
        val threeprime_whitelist
        val barcode_correction_method
        val barcode_correction_percentile

    output:
        tuple val(sample_id), path("${sample_id}.corrected.sorted.bam"), emit: barcode_corrected_tuple
        path "*"


    script:
        """
        # Correct step
        if [[ "${barcode_correction_method}" == "percentile" ]]; then
            isoseq correct --method "${barcode_correction_method}" --percentile "${barcode_correction_percentile}" --barcodes "${threeprime_whitelist}" "${refined_reads_bam}" "${sample_id}.corrected.bam"
        elif [[ "${barcode_correction_method}" == "knee" ]]; then
            isoseq correct --method "${barcode_correction_method}" --barcodes "${threeprime_whitelist}" "${refined_reads_bam}" "${sample_id}.corrected.bam"
        else
            echo "Invalid barcode correction method: ${barcode_correction_method}" >&2
            exit 1
        fi

        # Sort step
        samtools sort -t CB "${sample_id}.corrected.bam" -o "${sample_id}.corrected.sorted.bam"
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


process dedup_reads {
    label 'deduplication'

    input:
        tuple val(sample_id), path(barcode_corrected_bam)
        val barcode_correction_method
        val barcode_correction_percentile


    output:
        tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: dedup_tuple
        path "*"

    script:
    """
    isoseq groupdedup --keep-non-real-cells ${barcode_corrected_bam} ${sample_id}.dedup.bam
    if [[ "${barcode_correction_method}" == "percentile" ]]; then
      isoseq bcstats --method ${barcode_correction_method} --percentile ${barcode_correction_percentile} --json ${sample_id}.dedup.json -o ${sample_id}.dedup.tsv ${sample_id}.dedup.bam
    elif [[ "${barcode_correction_method}" == "knee" ]]; then
      isoseq bcstats --method ${barcode_correction_method} --json ${sample_id}.dedup.json -o ${sample_id}.dedup.tsv ${sample_id}.dedup.bam
    else
        echo "Invalid barcode correction method: ${barcode_correction_method}" >&2
        exit 1
    fi;
    samtools index ${sample_id}.dedup.bam

    """
}
