include { MTX_SUBSET } from '../../modules/mtx_subset.nf'
include {collect_mtx_as_h5ad as collect_filtered_isoforms_mtx_as_h5ad} from '../../modules/isoquant.nf'

workflow mtx_subset_wf {

  take:
    mtx_ch
    subset_f
  main:
    filtered_mtx_ch=MTX_SUBSET(mtx_ch,subset_f)
    filtered_h5ad_ch=collect_filtered_isoforms_mtx_as_h5ad(filtered_mtx_ch.collect(),"isoforms")




  emit:
    isoform_mtx=filtered_mtx_ch
    h5ad_file=filtered_h5ad_ch
}
