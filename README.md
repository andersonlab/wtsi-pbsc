# `wtsi-pbsc`: a Pacbio long-read single-cell RNA-seq processing pipeline for large-scale analyses based on Nextflow

This pipeline based on Nextflow performs common large-scale preprocessing for Pacbio single-cell RNA-seq data prepared using the Kinnex library preparation kit. Although most of the tools provided for the initial preprocessing steps are fairly easy to run, I wrote this in order to make the analysis of a large number of samples easier with nextflow.

## Specifying pipeline parameters
To run the pipeline you need to create:
 - A parameter file `params.yaml`: a YAML file that contains these required parameters to run the different steps. Not that each module will require different sets of parameters (detailed in each section below). These YAML syntax to specify these parameters:
```yaml
#Comment
param_path: "/path/to/file.txt"
param_integer: 20
```
 - An executor configuration file `exec.config`: specifies the runtime parameters such as memory and cpus to ask for for each step as well as the executor, queue name and error strategy. You can find more on the Nextflow documentation website [here](https://www.nextflow.io/docs/latest/config.html). We also provide a file `sanger.config` with the parameters that we used for isogut, although the appropriate runtime parameters may need to change per dataset and per execution environemnt.

## Components of `wtsi-pbsc`
`wtsi-pbsc` consists of four modules that represent different stages of pre-processing. Each module can be run indpendently provided the input files and parameters are correctly specified (e.g. using Nextflow option `entry fltnc`). Each step needs to be run with a set of parameters specified in the `parameters.yaml` file.
### 1- Step 1: HiFi reads to FLTNC reads
Using module `fltnc`. `fltnc` takes BAM files with HiFi reads (usually what you get from sequencing) and produces **f**ull-**l**ength **t**agged **n**on-**c**oncatemaer reads (FLTNC reads).
### 2- Step 2: Barcode correction and PCR deduplication
### 3- Step 3: Mapping to the genome using pbmm2 (which is a Pacbio wrapper around minimap2)
### 4- Step 4: Quantification using IsoQuant
As Quantification can be extremely time and memory-consuming, we implement several parallelisation strategies to speed up the process. 
#### No parallelisation: 
IsoQuant runs on all samples jointly and using all regions of the genome together. Suitable for a small number of samples (default IsoQuant behaviour)
#### By chromosome: 
IsoQuant runs on all samples jointly but splitting by chromosome. Suitable for a moderate number of samples (e.g. < 20)
#### By chunk: 
IsoQuant runs on all samples jointly but further splitting each chromosome into smaller chunks. Splitting points within regions of zero-coverage across all samples. Suitable for a large number of samples (e.g. < 60).
#### Two-pass approach: 
A two-pass approach where IsoQuant first quantifies known isoforms on a sample-by-sample basis (i.e. using IsoQuant's command-line argument `--no_model_construction`). Second, IsoQuant pools all remaining reads turning on isoform discovery (i.e. `inconsistent`,`inconsistent_ambiguous` and `inconsistent_non_intronic`).


### Parameters in `params.yaml`:

| Parameter                        | Description                                                                                                                                                             |
|-----------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `input_samples_path`              | Path to a comma-separated file with a header `sample_id`,`long_read_path`. Each entry should contain the `sample_id` and the path to the input BAM file with HiFi reads.|
| `exclude_samples`                 | File containing at least one column (`sample_id`) and listing the samples to be excluded from analysis.                                                                  |
| `skera_primers`                   | Path to primers that link concatenated reads. Required for `skera` to de-concatenate reads. Download [here](https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-MAS_adapters/MAS-Seq_Adapter_v1/mas16_primers.fasta). |
| `tenx_primers`                    | Path to 10X 3'/5' primers. Required for `lima` to remove 10X primers. Download [here](https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-10x_primers/). |
| `threeprime_whitelist`            | 10X cell barcode whitelist. Download [here](https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-10x_barcodes).                                                |
| `barcode_correction_method`       | Method used by `isoseq` to correct cell barcodes. Available options: `percentile`, `knee`. Default is `percentile`. Note that this method works better in our data.   |
| `barcode_correction_percentile`   | If `barcode_correction_method` is specified, use this to set a threshold for cell calling.                                                                            |
| `gtf_f`                           | GTF file used for quantification. This pipeline was tested with GENCODE v46.                                                                                           |
| `genome_fasta_f`                  | Genome in FASTA format.                                                                                                                                                 |
| `min_polya_length`                | Minimum polyA tail length. Any reads with fewer A in their polyA tails will be removed. Recommended: 20.                                                              |
| `results_output`                  | Path to output directory. Note that BAM and other types of files will be stored inside a subdirectory `qc`.                                                            |
| `chunks`                          | Number of chunks used for the parallelisation of IsoQuant.                                                                                                             |
| `isoquant_exclusion_regions_bed`  | Regions to exclude from quantification. In the data directory, V(D)J regions are provided as they are problematic for quantification in a large number of samples. Any BED file can be provided. |
















