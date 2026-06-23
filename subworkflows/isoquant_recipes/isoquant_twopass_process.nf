include {genedb_perChr_wf} from '../isoquant_components/genedb.nf'
include {isoquant_twopass_chunked_wf; isoquant_chrM;collect_output_wf} from './isoquant_twopass.nf'
include {chroms} from '../core/chroms.nf'

workflow ISOQUANT_TWOPASS_PROCESS {

  take:
    fullBam_ch
  main:
    if (fullBam_ch == 'independent workflow'){
      def samplesheet_header = file(params.input_samples_path).readLines()[0].split(',') as List
      if ('mapped_bam' in samplesheet_header) {
        Channel
          .fromPath(params.input_samples_path)
          .splitCsv(sep: ',', header: true)
          .map { it -> [it.sample_id, it.mapped_bam, it.mapped_bam + ".bai"] }
          .set { fullBam_ch }
      } else {
        Channel
          .fromPath(params.input_samples_path)
          .splitCsv(sep: ',', header: true)
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

        def deconv_dir_path = file("${params.results_output}deconvolution/bam/")
        if (deconv_dir_path.exists()) {
          Channel
            .fromPath(params.input_samples_path)
            .splitCsv(sep: ',', header: true)
            .filter { it -> it.nr_samples_multiplexed != '1' }
            .flatMap { it ->
              def sample_id = it.sample_id
              def nr_donors = it.nr_samples_multiplexed as int
              (0..<nr_donors).collect { n ->
                def bam_path = "${params.results_output}deconvolution/bam/${sample_id}_donor${n}.nosupplementary.bam"
                tuple("${sample_id}_donor${n}" as String, file(bam_path), file("${bam_path}.bai"))
              }
            }
            .set { split_bams }
          fullBam_ch = fullBam_ch.mix(split_bams)
        }
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
    preprocessed_bam_perChr_ch=chrom_ch.combine(fullBam_ch)

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
    (nochrM_output_chs.assignment_reads).concat(chrM_output_chs.assignment_reads),
    (nochrM_output_chs.transcriptmodel_reads).concat(chrM_output_chs.transcriptmodel_reads),
    (nochrM_output_chs.corrected_reads).concat(chrM_output_chs.corrected_reads)
    )
}