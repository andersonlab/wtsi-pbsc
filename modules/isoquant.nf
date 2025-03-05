process create_genedb_fasta_perChr {
  label 'micro_job'

  input:
      tuple val(chrom), path(gtf_f), path(fasta_f)
  output:
      tuple val(chrom), path("${chrom}.gtf"),path("${chrom}.gtf.db"), path("${chrom}.fa"), path("${chrom}.fa.fai")
  script:
  """
  awk -v chrom="${chrom}" '{if(\$1 ~ /^#/){print; next} else {if (\$1==chrom){print}} }' ${gtf_f} > "${chrom}.gtf"
  python ${baseDir}/scripts/create_genedb.py -g "${chrom}.gtf" -o "${chrom}.gtf.db"
  samtools faidx ${fasta_f} ${chrom} > "${chrom}.fa"
  samtools faidx "${chrom}.fa"
  """
}

process preprocess_bam {

  label 'micro_multithread_job'

  input:
      tuple val(sample_id), path(bam), path(bai)
  output:
      tuple val(sample_id), path("${sample_id}.mapped.realcells_only.processed.bam"), path("${sample_id}.mapped.realcells_only.processed.bam.bai")
  script:
  """
  sample_id="${sample_id}"

  #Removes supplementary alignments (note this may affect the detection of chimeric/fusion genes)
  #Appends sample name to CB tag
  samtools view -@ ${task.cpus} -h -F 2048 ${bam} |\\
  awk -v sample=\${sample_id} 'BEGIN {OFS="\t"}{if (\$1 ~ /^@/) {print; next}for(i=12; i<=NF; i++) {if (\$i ~ /^CB:Z:/) {\$i = \$i"_"sample}}print}' |\\
  samtools view -@ ${task.cpus} -h -bo "\${sample_id}.mapped.realcells_only.processed.bam" - ;
  samtools index -@ ${task.cpus} "\${sample_id}.mapped.realcells_only.processed.bam" ;
  """

}


process find_mapped_and_unmapped_regions_per_sampleChrom {
  label 'medium_job'

  input:
      tuple val(sample_id), val(chrom), path(bam), path(bai)
      path chrom_sizes_f
    output:
        tuple val(chrom), val(sample_id), path("${sample_id}_unmapped_regions.${chrom}.bed"), path("${sample_id}_mapped_regions.${chrom}.bed")
    script:
    """
    mapped_regions_f="${sample_id}_mapped_regions.${chrom}.bed"
    unmapped_regions_f="${sample_id}_unmapped_regions.${chrom}.bed"

    samtools view -F 4 -b ${bam} ${chrom} | bedtools bamtobed -i - | bedtools merge > \$mapped_regions_f
    grep -w $chrom ${chrom_sizes_f} | bedtools complement -i \$mapped_regions_f -g stdin > \$unmapped_regions_f
    """
}

process acrossSamples_mapped_unmapped_regions_perChr {
  label 'mini_job'

  input:
      tuple val(chrom), val(sample_ids), path(unmapped_beds), path(mapped_beds)
    output:
        tuple val(chrom), path("split_points_${chrom}.bed"), emit: unmapped_bed
        tuple val(chrom), path("mapped_regions_${chrom}.bed"), emit: mapped_bed
    script:
    """
    mapped_beds=(${mapped_beds.join(' ')})
    unmapped_beds=(${unmapped_beds.join(' ')})

    mapped_output_f=mapped_regions_"${chrom}".bed
    unmapped_output_f=split_points_"${chrom}".bed

    count_all_files="\${#mapped_beds[@]}";

    cat "\${mapped_beds[@]}" | sort -k1,1V -k2,2n | bedtools merge -i stdin > \$mapped_output_f;
    bedtools multiinter -i \${unmapped_beds[@]} | awk -v count_all_files="\${count_all_files}" '\$4==count_all_files' | cut -f1,2,3 > \$unmapped_output_f;
    """
}

process suggest_splits_binarySearch_v2 {
  label 'mini_job_local'

    input:
      tuple val(chrom), val(sample_ids), path(bams), path(bais), path(mapped_regions_bed)
      val chunks
      path chrom_sizes_f
    output:
      path "bams.txt"
      ///tuple val(chrom), path("suggested_splits_onebased_coords.${chrom}.bed"), path("suggested_splits_onebased_coords.${chrom}.list"), path("suggested_splits.${chrom}.bed")
    script:
    """
    printf "${bams.join('\n')}" > bams.txt
    """
}

process suggest_splits_binarySearch {
  label 'mini_job'

    input:
      tuple val(chrom), val(sample_ids), path(bams), path(bais), path(unmapped_regions_bed)
      val chunks
      path chrom_sizes_f
    output:
      tuple val(chrom), path("suggested_splits_onebased_coords.${chrom}.bed"), path("suggested_splits_onebased_coords.${chrom}.list"), path("suggested_splits.${chrom}.bed")
    script:
    """
    printf "${bams.join('\n')}" > bams.txt
    python ${baseDir}/dev/split_chr.py  -c ${chunks} -b bams.txt -r "${chrom}" -s ${unmapped_regions_bed} -z ${chrom_sizes_f} -o suggested_splits."${chrom}".bed
    awk -F "\t" '{print \$1"\t"(\$2+1)"\t"\$3}' suggested_splits."${chrom}".bed > suggested_splits_onebased_coords."${chrom}".bed
    awk -F "\t" '{print \$1":"\$2"-"\$3}' suggested_splits_onebased_coords."${chrom}".bed  > suggested_splits_onebased_coords."${chrom}".list
    rm bams.txt
    """
}

process split_bams {
  label 'micro_multithread_job'
  input:
    tuple val(chrom), val(sample_ids), path(bams), path(bais), val(formattedRegion), val(programmaticRegion)
  output:
    tuple val(chrom), val(sample_ids), path(outputBams), path(outputBais), val(formattedRegion), val(programmaticRegion),path("${programmaticRegion}_total_count.csv")
  script:

  outputBams = sample_ids.collect { sample_id -> sample_id+".${programmaticRegion}.mapped.realcells_only.processed.bam" }
  outputBais = sample_ids.collect { sample_id -> sample_id+".${programmaticRegion}.mapped.realcells_only.processed.bam.bai" }

  """
  #setting local bash variables
  bams=(${bams.join(' ')})
  sample_ids=(${sample_ids.join(' ')})
  programmaticRegion="${programmaticRegion}"
  formattedRegion="${formattedRegion}"
  chrom="${chrom}"

  total_count=0
  num_samples="\${#sample_ids[@]}"
  for i in \$(seq 0 \$((\$num_samples-1)) ); do
    sample_id="\${sample_ids[\$i]}"
    input_bam="\${bams[\$i]}"
    output_bam="\${sample_id}.\${programmaticRegion}.mapped.realcells_only.processed.bam"

    samtools view -@ ${task.cpus} -h "\${input_bam}" "\${formattedRegion}" |\
    awk -v chrom="\${chrom}" '{ if (\$1=="@SQ") {if(\$2=="SN:"chrom) {print \$0}} else{print \$0} }' |\
    samtools view -@ ${task.cpus} -h -bo "\${output_bam}" -
    samtools index -@ ${task.cpus} "\${output_bam}";
    count=\$(samtools view -@ ${task.cpus} -c \$output_bam)
    total_count=\$((\$total_count+\$count))
  done;
  echo "\${formattedRegion},\${total_count}" > "\${programmaticRegion}_total_count.csv"
  """


}


process run_isoquant_chunked {
    label 'isoquant_chunked'
    memory {
      def numReads = numReads.toInteger()
      def baseMemGB=3
      def additionalmemGB = numReads <= 100 ? 0 : numReads <= 1000 ? 2 : numReads <= 100000 ? 10 : numReads <= 1000000 ? 20 : numReads <= 10000000 ? 50 : numReads <= 100000000 ? 200 : 400
      return ((baseMemGB+additionalmemGB + (0.25 * additionalmemGB * (task.attempt-1) )).toInteger().toString()) + '.GB'
    }

    input:
        tuple val(chrom), val(sample_ids), path(bams), path(bais), val(formattedRegion), val(programmaticRegion), val(numReads), path(genedb), path(fasta), path(fai)




    output:
        tuple val(chrom), val(programmaticRegion), path("${programmaticRegion}/"), val("${numReads}"),path("custom_isoquantlog_${programmaticRegion}.out")


    script:
    """
    echo "${chrom},${formattedRegion},${programmaticRegion},${numReads},${task.memory}" > custom_isoquantlog_${programmaticRegion}.out


    isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${bams.join(' ')} --labels ${sample_ids.join(' ')} --data_type pacbio_ccs -o ${programmaticRegion} -p ${programmaticRegion} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format linear --bam_tags CB --no_secondary --clean_start


    """
}

process replace_novel_names {
    label 'micro_job'


    input:
      tuple val(chrom), val(programmaticRegion), path(isoquant_output)

    output:
        tuple val(chrom), val(programmaticRegion), path("${programmaticRegion}/${programmaticRegion}_renamed/")


    script:
    """
    input_dir=${isoquant_output}/${programmaticRegion}/
    output_dir=${isoquant_output}/${programmaticRegion}_renamed/
    mkdir -p \$output_dir

    transcriptgenefix_file_suffixes=(\\
    .transcript_model_counts.tsv \\
    .transcript_model_grouped_counts_linear.tsv \\
    .transcript_model_grouped_counts.tsv \\
    .transcript_model_grouped_tpm.tsv \\
    .transcript_model_reads.tsv.gz \\
    .transcript_models.gtf \\
    .transcript_model_tpm.tsv \\
    .transcript_tpm.tsv \\
    .transcript_counts.tsv \\
    .transcript_grouped_counts_linear.tsv \\
    .transcript_grouped_counts.tsv \\
    .transcript_grouped_tpm.tsv \\
    .novel_vs_known.SQANTI-like.tsv \\
    .extended_annotation.gtf \\
    .gene_counts.tsv \\
    .gene_grouped_counts_linear.tsv \\
    .gene_grouped_counts.tsv \\
    .gene_grouped_tpm.tsv \\
    .gene_tpm.tsv \\
    )

    asis_file_suffixes=(\\
    .intron_grouped_counts.tsv \\
    .intron_counts.tsv \\
    .exon_counts.tsv \\
    .exon_grouped_counts.tsv \\
    .corrected_reads.bed.gz \\
    .read_assignments.tsv.gz \\
    )

    exonfix_file_suffixes=(\\
    .extended_annotation.gtf \\
    .transcript_models.gtf \\
    )


    for suffix in \${transcriptgenefix_file_suffixes[@]}; do
      input_f="\${input_dir}${programmaticRegion}\${suffix}"
      output_f="\${output_dir}${programmaticRegion}\${suffix}"
      sed -E "s/(transcript[0-9]+)\\.([^.]+)\\.([^.]+)/\\1.${programmaticRegion}.\\3/g; s/(novel_gene)_([^_]+)_([0-9]+)/\\1_${programmaticRegion}_\\3/g" \${input_f} > \${output_f}
    done;

    for suffix in \${asis_file_suffixes[@]}; do
      input_f="\${input_dir}${programmaticRegion}\${suffix}"
      output_f="\${output_dir}${programmaticRegion}\${suffix}"
      cp \${input_f} \${output_f}
    done;

    #We also need to fix exon_ids in both extended_annotation and transcript_models GTFs
    for suffix in \${exonfix_file_suffixes[@]}; do
      #saving exonfixed GTFs as tmp file
      output_f_noexonfix="\${output_dir}${programmaticRegion}\${suffix}";
      output_f_withexonfixtmp="\${output_dir}${programmaticRegion}\${suffix}.tmp";
      bash ${baseDir}/scripts/fix_exon_ids.sh "\${output_f_noexonfix}" "\${output_f_withexonfixtmp}" "${programmaticRegion}";
      #reverting to original name
      rm \${output_f_noexonfix}
      mv \${output_f_withexonfixtmp} \${output_f_noexonfix}
    done;

    #This should later be added to isoquant_chunked process
    #There is a bug in Isoquant where if i include only inconsistent reads it still generates known isoforms in the transcript_model_linear_counts.tsv file
    output_suffix_withknown=.transcript_model_grouped_counts_linear.tsv
    output_suffix_noknown=.transcript_model_grouped_counts_linear.noknown.tsv
    output_f_withknown="\${output_dir}${programmaticRegion}\${output_suffix_withknown}"
    output_f_noknown="\${output_dir}${programmaticRegion}\${output_suffix_noknown}"
    grep -v "^ENST" \${output_f_withknown} > \${output_f_noknown}

    """
}

///////////////////////////////////////////////////
//////////IsoQuant output collection scripts////////
///////////////////////////////////////////////////
process collect_counts_as_mtx {
    label 'massive_job'

    input:
        path(isoquant_linear_count_files)

    output:
        path("barcodes.tsv")
        path("genes.tsv")
        path("matrix.mtx")

    script:
    """
    python ${baseDir}/scripts/convert_linear_counts_to_mtx.py -i ${isoquant_linear_count_files.join(' ')}
    """
}

process collect_counts_as_mtx_perChr {
    label 'big_job'

    input:
        tuple val(chrom), path(isoquant_linear_count_files)

    output:
        path("${chrom}"), emit: chrom_mtx
        path("${chrom}/barcodes.tsv")
        path("${chrom}/genes.tsv")
        path("${chrom}/matrix.mtx")

    script:
    """

    mkdir -p ${chrom}
    python ${baseDir}/scripts/convert_linear_counts_to_mtx.py -i ${isoquant_linear_count_files.join(' ')} -d ${chrom}/
    """
}

process collect_mtx_as_h5ad {
    label 'big_job'

    input:
        path(mtx_files)
        val(prefix)

    output:
        path("${prefix}.h5ad"), emit: h5ad_file
    script:
    """

    python ${baseDir}/scripts/mtx_to_hda5.py -i ${mtx_files.join(' ')} -p ${prefix}
    """
}

///Note that there shouldn't be any duplicates in GTFs for non-overlapping regions
process collect_gtfs {
    label 'big_job'

    input:
        path(query_gtf_files)
        path(ref_gtf_f)
        path(mtx_isoform_fs, stageAs: 'isoforms/isoforms?.tsv')


    output:
        path("extended_annotation.gtf")
        path("transcript_models.gtf")

    script:
    """
    query_gft_fs=(${query_gtf_files.join(' ')})
    mtx_isoform_fs=(${mtx_isoform_fs.join(' ')})

    for f in isoforms/isoforms*.tsv; do cut -f1 \$f; done | sort | uniq > all_features.csv
    for f in "\${query_gft_fs[@]}"; do echo \$f; done > query_gtf_files.txt

    python ${baseDir}/scripts/collect_gtfs.py -Q query_gtf_files.txt -r ${ref_gtf_f} -o extended_annotation.gtf
    echo "Finished collecting extended annotation GTF"
    python ${baseDir}/scripts/create_genedb.py -g extended_annotation.gtf -o extended_annotation.gtf.db
    echo "Finished creating extended annotation DB"
    python ${baseDir}/scripts/db_subset.py -d extended_annotation.gtf.db -i all_features.csv -o transcript_models.gtf
    echo "Finished subsetting DB"
    """
}

process format_intron_exon_grouped_counts {
    label 'small_job'

    input:
        tuple val(prefix), path(input_f)


    output:
        path("${prefix}.include_counts.tsv"), emit: include_counts
        path("${prefix}.exclude_counts.tsv"), emit: exclude_counts


    script:
    """
    python ${baseDir}/scripts/collect_intron_exon_grouped_counts.py -i ${input_f} -p ${prefix}
    """
}
process format_intron_exon_grouped_counts_perChr {
    label 'small_job'

    input:
        tuple val(chrom), val(prefix), path(input_f)
    output:
        tuple val(chrom), path("${prefix}.include_counts.tsv"), path("${prefix}.exclude_counts.tsv")

    script:
    """
    python ${baseDir}/scripts/collect_intron_exon_grouped_counts.py -i ${input_f} -p ${prefix}
    """
}




/////////Split by chromosome only not by chunk///////////
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
process run_isoquant_perChr {
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

/////////Two-pass IsoQuant///////////
process run_isoquant_firstPass {
label 'big_job'

  input:
      tuple val(chrom), val(sample_id), path(bam), path(bai), path(genedb), path(fasta), path(fai)
  output:
      tuple val(chrom), val(sample_id), path("${sample_id}/"),path("${sample_id}/${sample_id}.${chrom}/${sample_id}.${chrom}.read_assignments.tsv.gz"), path(bam)


  script:
  """
  isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${bam} --labels ${sample_id} --data_type pacbio_ccs -o ${sample_id} -p ${sample_id}.${chrom} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format linear --bam_tags CB --no_secondary --debug --no_model_construction
  """
}

process create_model_construction_bam {
label 'mini_job_local'

  input:
      tuple val(chrom), val(sample_id),path(read_assignment_f), path(bam)
  output:
      tuple val(chrom), val(sample_id), path("${sample_id}.${chrom}.mapped.realcells_only.processed.model_construction_reads.bam"), path("${sample_id}.${chrom}.mapped.realcells_only.processed.model_construction_reads.bam.bai")
  script:
  """
  model_construction_reads_list="${sample_id}.${chrom}.model_construction_reads.txt"
  model_construction_bam="${sample_id}.${chrom}.mapped.realcells_only.processed.model_construction_reads.bam"

  zcat ${read_assignment_f} | tail -n+4 | awk '{if( (\$6=="intergenic") || (\$6=="inconsistent_ambiguous") || (\$6=="inconsistent") || (\$6=="inconsistent_non_intronic"))print \$1}' | sort | uniq > \${model_construction_reads_list}

  samtools view -N \${model_construction_reads_list} -h -bo \${model_construction_bam} ${bam}
  samtools index \${model_construction_bam}
  """
}


process run_isoquant_secondPass {
label 'mini_job_local'

  input:
      tuple val(chrom), val(sample_ids), path(model_consutrciont_bams),path(model_consutrciont_bais), path(genedb), path(fasta), path(fai)
  output:
      tuple val(chrom), path("${chrom}/")
  script:
  """
    isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${model_consutrciont_bams.join(' ')} --labels ${sample_ids.join(' ')} --data_type pacbio_ccs -o ${chrom} -p ${chrom} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format linear --bam_tags CB --no_secondary --debug
  """
}
