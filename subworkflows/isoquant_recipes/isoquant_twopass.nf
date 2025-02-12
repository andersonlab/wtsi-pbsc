include {run_isoquant_firstPass; create_model_construction_bam; run_isoquant_chunked} from '../../modules/isoquant.nf'
include {find_mapped_and_unmapped_regions_per_sampleChrom; acrossSamples_mapped_unmapped_regions_perChr; suggest_splits_binarySearch; split_bams_perChunk} from '../../modules/smartSplit.nf'

workflow isoquant_twopass_perChr_wf {
  take:
    isoquant_preprocess_bam_perChr_ch
    chrom_genedb_fasta_chr_ch
  main:
    isoqunat_firsspass_input_ch=isoquant_preprocess_bam_perChr_ch
    .combine(chrom_genedb_fasta_chr_ch,by:0)

    isoquant_firstpass_output_ch=run_isoquant_firstPass(isoqunat_firsspass_input_ch)
    isoquant_firstpass_output_ch
    .map{ chrom,sample_id,isoquant_output_dir,read_assignment_f,bam -> [chrom,sample_id,read_assignment_f,bam] }
    .set{ model_construction_bam_input_ch }
    model_construction_bam_ch=create_model_construction_bam(model_construction_bam_input_ch)

    model_construction_bam_ch
    .groupTuple(by:0)
    .combine(chrom_genedb_fasta_chr_ch,by:0)
    .set{ isoquant_secondpass_input_ch }

    isoquant_secondpass_output_ch=run_isoquant_secondPass(isoquant_secondpass_input_ch)



  emit:
    isoquant_secondpass_output_ch
}



workflow isoquant_twopass_chunked_wf {
  take:
    isoquant_preprocess_bam_perChr_ch
    chrom_genedb_fasta_chr_ch
    chrom_sizes_f
    chunks
  main:
    isoqunat_firsspass_input_ch=isoquant_preprocess_bam_perChr_ch
    .combine(chrom_genedb_fasta_chr_ch,by:0)

    isoquant_firstpass_output_ch=run_isoquant_firstPass(isoqunat_firsspass_input_ch)
    isoquant_firstpass_output_ch
    .map{ chrom,sample_id,isoquant_output_dir,read_assignment_f,bam -> [chrom,sample_id,read_assignment_f,bam] }
    .set{ model_construction_bam_input_ch }
    model_construction_bam_ch=create_model_construction_bam(model_construction_bam_input_ch)

    ///sharding
    mapped_unmapped_regions_tuple_ch=find_mapped_and_unmapped_regions_per_sampleChrom(model_construction_bam_ch,chrom_sizes_f)
    mapped_unmapped_regions_tuple_ch.groupTuple().set{mapped_unmapped_regions_groupedTuple_ch}
    acrossSamples_mapped_unmapped_regions_perChr_ch=acrossSamples_mapped_unmapped_regions_perChr(mapped_unmapped_regions_groupedTuple_ch)


    /// Merging grouped BAM tuples with mapped/unmapped BED files
    model_construction_bam_ch.groupTuple(by:0).combine(acrossSamples_mapped_unmapped_regions_perChr_ch.mapped_bed,by:0).set { modelconstructionBam_mappedbed_tuples_ch }
    model_construction_bam_ch.groupTuple(by:0).combine(acrossSamples_mapped_unmapped_regions_perChr_ch.unmapped_bed,by:0).set { modelconstructionBam_unmappedbed_tuples_ch }

    suggested_splits_ch=suggest_splits_binarySearch(modelconstructionBam_unmappedbed_tuples_ch,chunks,chrom_sizes_f)


    suggested_splits_ch.flatMap { tuple ->
        def chrom = tuple[0]     // Extract chrom from the tuple
        def filePath = tuple[3]  // Get the last item in the tuple (path to file)
        file(filePath).text.split('\n') // Read file content and split into lines
            .findAll { it }             // Remove empty lines
            .collect { line ->          // Format each line and pair it with chrom
                def cols = line.split('\t') // Split line into columns by tab
                def formattedRegion = cols[0]+":"+cols[1]+"-"+cols[2] // Region formatting
                def programmaticRegion=cols[0]+"_"+cols[1]+"_"+cols[2]
                [chrom, formattedRegion,programmaticRegion] // Return chrom and formatted line as a pair
            }
    }.set {chrom_region_ch}

    /// Splitting BAMs according to suggested regions
    model_construction_bam_ch
      .groupTuple(by:0)
      .combine(chrom_region_ch,by:0)
      .set {modelconstructionBam_Region_groupedTuple_ch}
    modelconstructionRegionBam=split_bams_perChunk(modelconstructionBam_Region_groupedTuple_ch,"mapped.realcells_only.processed.model_construction_reads")


    ///Runnning second pass
    modelconstructionRegionBam
    .map {chrom, sample_ids, bams, bais, formattedRegion, programmaticRegion, counts_f  -> [chrom, sample_ids, bams, bais, formattedRegion, programmaticRegion, file(counts_f).splitCsv(header:false)[0][1] ] }
    .combine(chrom_genedb_fasta_chr_ch,by:0)
    .set {isoquant_chunked_input}



    isoquant_secondpass_output_ch=run_isoquant_chunked(isoquant_chunked_input)





  emit:
    isoquant_chunked_input
}
