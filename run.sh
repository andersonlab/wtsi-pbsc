#BSUB -J mo11_testrun
#BSUB -R "select[mem>5000] rusage[mem=5000]"
#BSUB -M 5000
#BSUB -n 2
#BSUB -a "memlimit=True"
#BSUB -o mo11_testrun_out.log
#BSUB -e mo11_testrun_err.log
#BSUB -q oversubscribed

#module load HGI/softpack/groups/team152/pblr-rnaseq/8
module load HGI/pipelines/yascp/1.9
module load HGI/softpack/groups/team152/pblr-rnaseq/11
#module load HGI/softpack/groups/hgi/pblr-rnaseq-test/1
module load HGI/common/nextflow/24.10.4
export NXF_OPTS="-Xms8G -Xmx16G"
exec_time=$(date -d "today" +"%Y%m%d-%H%M%S")
# workflow=fltnc
# workflow=correct_barcodes
workflow=map_pbmm2
# workflow=isoquant_twopass
#workflow=sqanti3
#workflow=isoquant_chunked_old
# workflow=full
workdir=$PWD/work
#workdir=/lustre/scratch126/humgen/teams/hgi/users/vo3/isoseq_test/isoseq_workdir_dedup_test

SINGULARITY_TMPDIR=$workdir/tmp
TEMP=$workdir/tmp
TMP_DIR=$workdir/tmp
export $SINGULARITY_TMPDIR
export $TEMP
export $TMP_DIR
nextflow run /lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/mo11/mo11_tmp_work/isoquant_tests/v2/merged/wtsi-pbsc/isoseq2.nf --input_samples_path /lustre/scratch127/humgen/teams_v2/hgi/mo11/tmp_projects/isoquant/input_test.csv -c input.nf -profile sanger -resume ${workdir} -entry ${workflow} -with-trace  -with-report report-${exec_time}.html -with-dag dag-${exec_time}.html > nextflow.nohup.log 2>&1 & 

