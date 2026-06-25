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
    // memory {
    //   def numReads = numReads.toInteger()
    //   def baseMemGB=3
    //   def additionalmemGB = numReads <= 100 ? 0 : numReads <= 1000 ? 2 : numReads <= 100000 ? 10 : numReads <= 1000000 ? 20 : numReads <= 10000000 ? 50 : numReads <= 100000000 ? 200 : 400
    //   return ((baseMemGB+additionalmemGB + (0.25 * additionalmemGB * (task.attempt-1) )).toInteger().toString()) + '.GB'
    // }

    input:
        tuple val(chrom), val(sample_ids), path(bams), path(bais), val(formattedRegion), val(programmaticRegion), path(genedb), path(fasta), path(fai)




    output:
        tuple val(chrom), val(programmaticRegion), path("${programmaticRegion}/")


    script:
    """
    isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${bams.join(' ')} --labels ${sample_ids.join(' ')} --data_type pacbio_ccs -o ${programmaticRegion} -p ${programmaticRegion} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format mtx --bam_tags CB --no_secondary --clean_start --polya_trimmed all --process_only_chr ${chrom}
    rm -f ${programmaticRegion}/${programmaticRegion}/${programmaticRegion}.extended_annotation.gtf
    """
}

process replace_novel_names {
    label 'micro_job'
    tag "${programmaticRegion}"

    input:
      tuple val(chrom), val(programmaticRegion), path(isoquant_tar)

    output:
      tuple val(chrom), val(programmaticRegion), path("${programmaticRegion}_renamed.tar")


    script:
    """
    tar -xzf ${isoquant_tar}
    input_dir=${programmaticRegion}/${programmaticRegion}/
    output_dir=${programmaticRegion}_renamed/
    mkdir -p \$output_dir

    # Load suffix lists from data files
    transcriptgenefix_suffixes=()
    while IFS= read -r line || [[ -n "\$line" ]]; do
      transcriptgenefix_suffixes+=("\$line")
    done < ${baseDir}/data/isoquant_transcriptgenefix_suffixes.txt
    exonfix_suffixes=()
    while IFS= read -r line || [[ -n "\$line" ]]; do
      exonfix_suffixes+=("\$line")
    done < ${baseDir}/data/isoquant_exonfix_suffixes.txt

    # Build lookup set for O(1) membership test
    declare -A transcriptgenefix_set
    for s in "\${transcriptgenefix_suffixes[@]}"; do
      transcriptgenefix_set["\$s"]=1
    done

    # Detect whether any novel transcript models exist in this chunk
    input_transcript_counts="\${input_dir}${programmaticRegion}.discovered_transcript_counts.tsv"
    if [[ -e "\$input_transcript_counts" ]] && [[ \$(wc -l < "\$input_transcript_counts") -gt 1 ]]; then
      HAS_NOVEL=true
    else
      HAS_NOVEL=false
    fi

    # Route every file: transcriptgenefix files get sed (or symlink if no novel),
    # everything else gets a symlink unconditionally
    for f in "\${input_dir}"*; do
      [[ -f "\$f" ]] || continue
      fname=\$(basename "\$f")
      suffix="\${fname#${programmaticRegion}}"
      output_f="\${output_dir}\${fname}"
      if [[ -n "\${transcriptgenefix_set["\$suffix"]+_}" ]]; then
        if [[ "\$HAS_NOVEL" == "true" ]]; then
          zcat -f "\$f" | sed -E "s/(transcript[0-9]+)\\.([^.]+)\\.([^.]+)/\\1.${programmaticRegion}.\\3/g; s/(novel_gene)_([^_]+)_([0-9]+)/\\1_${programmaticRegion}_\\3/g" > "\$output_f"
        else
          ln -s \$(realpath "\$f") "\$output_f"
        fi
      else
        ln -s \$(realpath "\$f") "\$output_f"
      fi
    done

    # Fix novel exon IDs in GTFs (only needed when novel models exist)
    if [[ "\$HAS_NOVEL" == "true" ]]; then
      for suffix in "\${exonfix_suffixes[@]}"; do
        output_f="\${output_dir}${programmaticRegion}\${suffix}"
        if [[ -e "\$output_f" ]]; then
          output_f_tmp="\${output_f}.tmp"
          bash ${baseDir}/scripts/fix_exon_ids.sh "\$output_f" "\$output_f_tmp" "${programmaticRegion}"
          rm "\$output_f"
          mv "\$output_f_tmp" "\$output_f"
        fi
      done
    fi

    # Generate .noknown.tsv files (filter known transcripts introduced by IsoQuant bug, fixed in 3.7.0)
    output_f_withknown="\${output_dir}${programmaticRegion}.discovered_transcript_counts.tsv"
    output_f_noknown="\${output_dir}${programmaticRegion}.discovered_transcript_counts.noknown.tsv"
    if [[ -e "\${output_f_withknown}" ]]; then
      grep -v -e "^ENST" -e "__ambiguous" -e "__no_feature" -e "__not_aligned" "\${output_f_withknown}" > "\${output_f_noknown}"
    fi

    output_f_withknown="\${output_dir}${programmaticRegion}.discovered_transcript_grouped_tag_CB_counts.linear.tsv"
    output_f_noknown="\${output_dir}${programmaticRegion}.discovered_transcript_grouped_tag_CB_counts.linear.noknown.tsv"
    if [[ -e "\${output_f_withknown}" ]]; then
      grep -v -e "^ENST" -e "__ambiguous" -e "__no_feature" -e "__not_aligned" "\${output_f_withknown}" > "\${output_f_noknown}"
    fi

    # --dereference resolves symlinks created for unchanged files so the tar contains real content
    tar --dereference -czf ${programmaticRegion}_renamed.tar ${programmaticRegion}_renamed/
    rm -rf ${programmaticRegion}_renamed/ ${programmaticRegion}/
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
    label 'counts_collect'
    tag "${chrom}"
    publishDir "${publish_dir}", mode: 'copy', overwrite: true

    input:
        // tars: one tar per sample/chunk; patterns: one wildcard per tar matching the count file inside
        // needs_chr_filter: parallel boolean list — true for per-sample whole-genome tars (filter rows
        //   to this chromosome using read_assignments), false for already-chromosome-scoped tars
        tuple val(chrom), path(tars), val(patterns), val(needs_chr_filter)
        val(publish_dir)

    output:
        path("${chrom}"), emit: chrom_mtx
        path("${chrom}/barcodes.tsv")
        path("${chrom}/genes.tsv")
        path("${chrom}/matrix.mtx")

    script:
    """
    tars=(${tars.join(' ')})
    patterns=(${patterns.join(' ')})
    needs_filter=(${needs_chr_filter.join(' ')})
    count_files=()

    for i in \$(seq 0 \$((\${#tars[@]}-1))); do
      out="count_\${i}.tsv"
      if [ "\${needs_filter[\$i]}" = "true" ]; then
        chr_features="chr_features_\${i}.txt"
        tar -xzf "\${tars[\$i]}" -O --wildcards "*.read_assignments.tsv.gz" | zcat | \
          awk -v chr="${chrom}" '
            BEGIN { FS="\\t" }
            /^#/ { next }
            !header_parsed {
              for (i=1;i<=NF;i++) { if (\$i=="chr") chr_col=i; if (\$i=="isoform_id") id_col=i; if (\$i=="gene_id") gene_col=i }
              header_parsed=1; next
            }
            chr_col>0 && \$chr_col==chr {
              if (id_col>0 && \$id_col!="" && \$id_col!=".") print \$id_col
              if (gene_col>0 && \$gene_col!="" && \$gene_col!=".") print \$gene_col
            }
          ' | sort | uniq > "\${chr_features}"
        tar -xzf "\${tars[\$i]}" -O --wildcards "\${patterns[\$i]}" > "\${out}.raw"
        awk 'NR==FNR{a[\$1]; next} FNR==1 || (\$1 in a)' "\${chr_features}" "\${out}.raw" > "\${out}"
        rm "\${chr_features}" "\${out}.raw"
      else
        tar -xzf "\${tars[\$i]}" -O --wildcards "\${patterns[\$i]}" > "\${out}"
      fi
      count_files+=("\${out}")
    done

    mkdir -p ${chrom}
    python ${baseDir}/scripts/convert_linear_counts_to_mtx.py -i "\${count_files[@]}" -d ${chrom}/
    printf '%s\0' "\${count_files[@]}" | xargs -0 rm -f
    """

}

process collect_mtx_as_h5ad {
    label 'counts_collect'
    tag "${prefix}"
    publishDir "${publish_dir}", mode: 'copy', overwrite: true

    input:
        path(mtx_files)
        val(prefix)
        val(publish_dir)

    output:
        path("${prefix}.h5ad"), emit: h5ad_file
    script:
    """
    python ${baseDir}/scripts/mtx_to_hda5.py -i ${mtx_files.join(' ')} -p ${prefix}
    """

}

///Note that there shouldn't be any duplicates in GTFs for non-overlapping regions
process collect_gtfs {
    label 'collect_gtfs'
    publishDir "${publish_dir}", mode: 'copy', overwrite: true
    input:
        path(renamed_tars)
        path(ref_gtf_f)
        path(mtx_isoform_fs, stageAs: 'isoforms/isoforms?.tsv')
        val(publish_dir)


    output:
        path("extended_annotation.gtf")
        path("transcript_models.gtf")
        path("extended_annotation.gtf.db")
        path("transcript_models.gtf.db")

    script:
    """
    renamed_tars=(${renamed_tars.join(' ')})

    for f in isoforms/isoforms*.tsv; do cut -f1 \$f; done | sort | uniq > all_features.csv

    gtf_files=()
    for i in \$(seq 0 \$((\${#renamed_tars[@]}-1))); do
      out="gtf_\${i}.gtf"
      tar -xzf "\${renamed_tars[\$i]}" -O --wildcards "*.transcript_models.gtf" > "\${out}"
      gtf_files+=("\${out}")
    done

    for f in "\${gtf_files[@]}"; do echo \$f; done > query_gtf_files.txt

    python ${baseDir}/scripts/collect_gtfs.py -Q query_gtf_files.txt -r ${ref_gtf_f} -o extended_annotation.gtf
    echo "Finished collecting extended annotation GTF"
    python ${baseDir}/scripts/create_genedb.py -g extended_annotation.gtf -o extended_annotation.gtf.db
    echo "Finished creating extended annotation DB"
    python ${baseDir}/scripts/db_subset.py -d extended_annotation.gtf.db -i all_features.csv -o transcript_models.gtf
    echo "Finished subsetting DB as GTF"
    python ${baseDir}/scripts/create_genedb.py -g transcript_models.gtf -o transcript_models.gtf.db
    echo "Finished converting GTF to DB"

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
    tag "${prefix}"

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
    isoquant.py --reference ${genome_fasta_f} --genedb ${gtf_f} --complete_genedb --sqanti_output --bam ${bams.join(' ')} --labels ${sample_ids.join(' ')} --data_type pacbio_ccs -o ${chrom} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format mtx --bam_tags CB --no_secondary --clean_start --polya_trimmed all --process_only_chr ${chrom}
    """
}

/////////Two-pass IsoQuant///////////
process run_isoquant_firstPass {
label 'isoquant_firstPass'
tag "${sample_id}__${chrom}"

  input:
      tuple val(chrom), val(sample_id), path(bam), path(bai)
      val(genedb)
      val(fasta)
      val(fai)
  output:
      tuple val(chrom), val(sample_id), path("${sample_id}.tar"), path(bam)

  script:
  """
  isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${bam} --labels ${sample_id} --data_type pacbio_ccs -o ${sample_id} -p ${sample_id}.${chrom} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format mtx --bam_tags CB --no_secondary --no_model_construction --polya_trimmed all --process_only_chr ${chrom}
  rm -f ${sample_id}/${sample_id}.${chrom}/${sample_id}.${chrom}.extended_annotation.gtf
  tar -czf ${sample_id}.tar ${sample_id}/
  rm -rf ${sample_id}/
  """
}
process run_isoquant_firstPass_perSample {
label 'isoquant_firstPass_perSample'
tag "${sample_id}"

  input:
      tuple val(sample_id), path(bam), path(bai)
      val(genedb)
      val(fasta)
      val(fai)
  output:
      tuple val(sample_id), path("${sample_id}.tar"), path(bam)

  script:
  """
  isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${bam} --labels ${sample_id} --data_type pacbio_ccs -o ${sample_id} -p ${sample_id} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format mtx --bam_tags CB --no_secondary --no_model_construction --polya_trimmed all
  rm -f ${sample_id}/${sample_id}/${sample_id}.extended_annotation.gtf
  tar -czf ${sample_id}.tar ${sample_id}/
  rm -rf ${sample_id}/
  """
}

////////////////////
///chrM processes///
////////////////////
process run_isoquant_firstPass_withmodelconstruction {
label 'isoquant_firstPass_withmodelconstruction'
tag "${sample_id}__${chrom}"

  input:
      tuple val(chrom), val(sample_id), path(bam), path(bai)
      val(genedb)
      val(fasta)
      val(fai)
  output:
      tuple val(chrom), val(sample_id), path("${sample_id}.tar"), path(bam)
  script:
  """
  isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${bam} --labels ${sample_id} --data_type pacbio_ccs -o ${sample_id} -p ${sample_id}.${chrom} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format mtx --bam_tags CB --no_secondary --polya_trimmed all --process_only_chr ${chrom}
  rm -f ${sample_id}/${sample_id}.${chrom}/${sample_id}.${chrom}.extended_annotation.gtf
  tar -czf ${sample_id}.tar ${sample_id}/
  rm -rf ${sample_id}/
  """
}

process replace_novel_names_firsPass_singlenovelname {
    label 'micro_job'
    tag "${sample_id}__${chrom}"

    input:
      tuple val(chrom), val(sample_id), path(isoquant_tar)

    output:
        tuple val(chrom), val(sample_id), path("${sample_id}.${chrom}_renamed.tar")


    script:
    """
    tar -xzf ${isoquant_tar}
    input_dir=${sample_id}/${sample_id}.${chrom}/
    output_dir=${sample_id}.${chrom}_renamed/
    mkdir -p \$output_dir

    # Load suffix lists from data files
    transcriptgenefix_suffixes=()
    while IFS= read -r line || [[ -n "\$line" ]]; do
      transcriptgenefix_suffixes+=("\$line")
    done < ${baseDir}/data/isoquant_transcriptgenefix_suffixes.txt
    exonfix_suffixes=()
    while IFS= read -r line || [[ -n "\$line" ]]; do
      exonfix_suffixes+=("\$line")
    done < ${baseDir}/data/isoquant_exonfix_suffixes.txt

    # Build lookup set for O(1) membership test
    declare -A transcriptgenefix_set
    for s in "\${transcriptgenefix_suffixes[@]}"; do
      transcriptgenefix_set["\$s"]=1
    done

    # Detect whether any novel transcript models exist in this chunk
    input_transcript_counts="\${input_dir}${sample_id}.${chrom}.discovered_transcript_counts.tsv"
    if [[ -e "\$input_transcript_counts" ]] && [[ \$(wc -l < "\$input_transcript_counts") -gt 1 ]]; then
      HAS_NOVEL=true
    else
      HAS_NOVEL=false
    fi

    # Route every file: transcriptgenefix files get sed (or symlink if no novel),
    # everything else gets a symlink unconditionally
    for f in "\${input_dir}"*; do
      [[ -f "\$f" ]] || continue
      fname=\$(basename "\$f")
      suffix="\${fname#${sample_id}.${chrom}}"
      output_f="\${output_dir}\${fname}"
      if [[ -n "\${transcriptgenefix_set["\$suffix"]+_}" ]]; then
        if [[ "\$HAS_NOVEL" == "true" ]]; then
          zcat -f "\$f" | sed -E "s/(transcript[0-9]+)\\.([^.]+)\\.([^.]+)/\\1.\\2_${sample_id}.\\3/g; s/(novel_gene)_([^_]+)_([0-9]+)/\\1_${sample_id}_\\3/g" > "\$output_f"
        else
          ln -s \$(realpath "\$f") "\$output_f"
        fi
      else
        ln -s \$(realpath "\$f") "\$output_f"
      fi
    done

    # Fix novel exon IDs in GTFs (only needed when novel models exist)
    if [[ "\$HAS_NOVEL" == "true" ]]; then
      for suffix in "\${exonfix_suffixes[@]}"; do
        output_f="\${output_dir}${sample_id}.${chrom}\${suffix}"
        if [[ -e "\$output_f" ]]; then
          output_f_tmp="\${output_f}.tmp"
          bash ${baseDir}/scripts/fix_exon_ids.sh "\$output_f" "\$output_f_tmp" "${sample_id}"
          rm "\$output_f"
          mv "\$output_f_tmp" "\$output_f"
        fi
      done
    fi

    # Generate .noknown.tsv (filter special-token rows introduced by IsoQuant bug, fixed in 3.7.0)
    output_f_withknown="\${output_dir}${sample_id}.${chrom}.discovered_transcript_counts.tsv"
    output_f_noknown="\${output_dir}${sample_id}.${chrom}.discovered_transcript_counts.noknown.tsv"
    if [[ -e "\${output_f_withknown}" ]]; then
      grep -v -e "__ambiguous" -e "__no_feature" -e "__not_aligned" "\${output_f_withknown}" > "\${output_f_noknown}"
    fi

    # --dereference resolves symlinks created for unchanged files so the tar contains real content
    tar --dereference -czf ${sample_id}.${chrom}_renamed.tar ${sample_id}.${chrom}_renamed/
    rm -rf ${sample_id}.${chrom}_renamed/ ${sample_id}/
    """
}
/////////////////////////////
//////END: chrM processes////
/////////////////////////////


process create_model_construction_bam {
label 'mini_job'

  input:
      tuple val(chrom), val(sample_id),path(read_assignment_f), path(bam)
  output:
      tuple val(chrom), val(sample_id), path("${sample_id}.${chrom}.model_construction_reads.bam"), path("${sample_id}.${chrom}.model_construction_reads.bam.bai")
  script:
  """
  model_construction_reads_list="${sample_id}.${chrom}.model_construction_reads.txt"
  model_construction_bam="${sample_id}.${chrom}.model_construction_reads.bam"

  zcat ${read_assignment_f} | tail -n+4 | awk '{if( (\$6=="intergenic") || (\$6=="inconsistent_ambiguous") || (\$6=="inconsistent") || (\$6=="inconsistent_non_intronic"))print \$1}' | sort | uniq > \${model_construction_reads_list}

  samtools view -N \${model_construction_reads_list} -h -bo \${model_construction_bam} ${bam}
  samtools index \${model_construction_bam}
  """
}


process create_model_construction_bam_perChr {
label 'medium_job'
tag "${chrom}"

  input:
      tuple val(chrom), val(sample_ids), path(firstpass_tars), path(bams)
  output:
      tuple val(chrom), path("${chrom}.model_construction_reads.bam"), path("${chrom}.model_construction_reads.bam.bai")
  script:
  """
  sample_ids=(${sample_ids.join(' ')})
  firstpass_tars=(${firstpass_tars.join(' ')})
  bams=(${bams.join(' ')})

  temp_bams=()
  for i in \$(seq 0 \$((\${#sample_ids[@]}-1))); do
    sample_id="\${sample_ids[\$i]}"
    tar_f="\${firstpass_tars[\$i]}"
    bam="\${bams[\$i]}"

    reads_list="\${sample_id}.${chrom}.model_construction_reads.txt"
    temp_bam="\${sample_id}.${chrom}.model_construction_reads.tmp.bam"

    tar -xzf "\${tar_f}" -O --wildcards "*.read_assignments.tsv.gz" | zcat | \
      awk -v chr="${chrom}" '
        BEGIN { FS="\\t" }
        /^#/ { next }
        !header_parsed {
          for (i=1;i<=NF;i++) { if (\$i=="chr") chr_col=i; if (\$i=="assignment_type") type_col=i }
          header_parsed=1; next
        }
        chr_col>0 && type_col>0 && \$chr_col==chr && (\$type_col=="intergenic"||\$type_col=="inconsistent_ambiguous"||\$type_col=="inconsistent"||\$type_col=="inconsistent_non_intronic") { print \$1 }
      ' | sort | uniq > "\${reads_list}"
    samtools view -N "\${reads_list}" -h -bo "\${temp_bam}" "\${bam}"
    temp_bams+=("\${temp_bam}")
  done

  samtools merge -f "${chrom}.model_construction_reads.bam" "\${temp_bams[@]}"
  samtools index "${chrom}.model_construction_reads.bam"
  printf '%s\0' "\${temp_bams[@]}" "\${firstpass_tars[@]}" | xargs -0 rm -f
  """
}


process run_isoquant_chunked_merged {
    label 'isoquant_chunked'
    tag "${programmaticRegion}"

    input:
        tuple val(chrom), path(bam), path(bai), val(formattedRegion), val(programmaticRegion)
        val(genedb)
        val(fasta)
        val(fai)

    output:
        tuple val(chrom), val(programmaticRegion), path("${programmaticRegion}.tar")

    script:
    """
    isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${bam} --labels ${programmaticRegion} --data_type pacbio_ccs -o ${programmaticRegion} -p ${programmaticRegion} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format mtx --bam_tags CB --no_secondary --clean_start --polya_trimmed all --process_only_chr ${chrom}
    rm -f ${programmaticRegion}/${programmaticRegion}/${programmaticRegion}.extended_annotation.gtf
    tar -czf ${programmaticRegion}.tar ${programmaticRegion}/
    rm -rf ${programmaticRegion}/
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
    isoquant.py --reference ${fasta} --genedb ${genedb} --complete_genedb --sqanti_output --bam ${model_consutrciont_bams.join(' ')} --labels ${sample_ids.join(' ')} --data_type pacbio_ccs -o ${chrom} -p ${chrom} --count_exons --check_canonical  --read_group tag:CB -t ${task.cpus} --counts_format mtx --bam_tags CB --no_secondary --polya_trimmed all --process_only_chr ${chrom}
  """
}
