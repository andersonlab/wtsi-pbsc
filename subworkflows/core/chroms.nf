workflow chroms {
  take:
    chromosomes_list

  main:
    ///Channel emitting chromosomes
    ///Testing
    /// def chromosomes_list = ['chr22']
    ///
    Channel.from(chromosomes_list).set {chrom_ch}

  emit:
    chrom_ch
}
