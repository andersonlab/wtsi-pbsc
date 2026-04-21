include {genedb_perChr_wf} from '../isoquant_components/genedb.nf'
include {preprocess_bam_perChr_wf} from '../isoquant_components/preprocess_bam.nf'
include {isoquant_twopass_chunked_wf; isoquant_chrM;collect_output_wf} from './isoquant_twopass.nf'
include {chroms} from '../core/chroms.nf'

workflow ISOQUANT_TWOPASS_PROCESS {

  take:
    fullBam_ch
  main:
    if (fullBam_ch == 'independent workflow'){
      Channel
        .fromPath(params.input_samples_path)
        .splitCsv(sep: ',', header: true)
        //.filter { it -> !(it.sample_id in excluded_samples_list) }
        .map { it ->
          def sample_id = it.sample_id
          def bam_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam"
          def bai_path = "${params.results_output}qc/mapped/${sample_id}.mapped.realcells_only.bam.bai"
          [sample_id, bam_path, bai_path]
        }
        .set { fullBam_ch }
      Channel
        .fromPath(params.input_samples_path)
        .splitCsv(sep: ',', header: true)
        .map { it -> [it.sample_id, it.nr_samples_multiplexed] } // Create a tuple with bam_path and sample_id
        .set { sampleNames_nrDons }

      fullBam_ch =  fullBam_ch.combine(sampleNames_nrDons, by: 0)
      fullBam_ch=fullBam_ch
      .filter { experiment, bam, bai, npooled -> npooled == '1' }
      .map { experiment, bam, bai, npooled -> [experiment, bam, bai] }

      def folder = file("${params.results_output}deconvolution/bam/")
      if (folder.exists()) {
        Channel
          .fromPath("${params.results_output}deconvolution/bam/*.bam")
          .map { bam_file -> 
              // Removes only the LAST _donorN suffix
              def name = bam_file.baseName.replaceAll(/_donor\d+$/, '')
              tuple(name, bam_file, file(bam_file.toString() + ".bai"))
          }
          .set { split_bams }
        split_bams =  split_bams.combine(sampleNames_nrDons, by: 0)
        split_bams = split_bams
        .filter { experiment, bam, bai, npooled -> npooled != '1' }
        .map { experiment, bam, bai, npooled -> [bam.simpleName, bam, bai] }
        fullBam_ch=fullBam_ch.mix(split_bams)
      }
    }
    def chromosomes_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                            'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                            'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY','chrM']

    if(params.chunks) {
    chrom_sizes_f=params.nf_basedir+"data/hg38.chrom.sizes"
    }

    chrom_ch=chroms(chromosomes_list)

    chrom_genedb_fasta_chr_ch=genedb_perChr_wf(chrom_ch,params.gtf_f,params.genome_fasta_f)
    preprocessed_bam_perChr_ch=preprocess_bam_perChr_wf(chrom_ch,fullBam_ch)

    //Processing chrM separately
    preprocessed_bam_perChr_ch.filter{chrom,sample_id,bam,bai -> chrom=="chrM"}.set{preprocessed_bam_chrM_ch}
    chrom_genedb_fasta_chr_ch.filter{chrom,gene_db,fa,fai -> chrom=="chrM"}.set{chrom_genedb_fasta_chrM_ch}
    chrM_output_chs=isoquant_chrM(preprocessed_bam_chrM_ch,chrom_genedb_fasta_chrM_ch,chrom_sizes_f)

    //Processing other chromosomes separately
    preprocessed_bam_perChr_ch.filter{chrom,sample_id,bam,bai -> chrom!="chrM"}.set{preprocessed_bam_nochrM_ch}
    chrom_genedb_fasta_chr_ch.filter{chrom,gene_db,fa,fai -> chrom!="chrM"}.set{chrom_genedb_fasta_nochrM_ch}
    nochrM_output_chs=isoquant_twopass_chunked_wf(preprocessed_bam_nochrM_ch,chrom_genedb_fasta_nochrM_ch,chrom_sizes_f,params.chunks)

    //Collecting output: MTX, GTF, reads
    collect_output_wf(
    (nochrM_output_chs.isoform_counts).concat(chrM_output_chs.isoform_counts),
    (nochrM_output_chs.gene_counts).concat(chrM_output_chs.gene_counts),
    (nochrM_output_chs.existing_gtf).concat(chrM_output_chs.existing_gtf),
    (nochrM_output_chs.extended_gtf).concat(chrM_output_chs.extended_gtf),
    (nochrM_output_chs.assignment_reads).concat(chrM_output_chs.assignment_reads),
    (nochrM_output_chs.transcriptmodel_reads).concat(chrM_output_chs.transcriptmodel_reads),
    (nochrM_output_chs.corrected_reads).concat(chrM_output_chs.corrected_reads)
    )
}