export NXF_OPTS="-Xms8G -Xmx16G"

#######################
######################
# Example run on LSF#
#Work directory for temporary files
workdir=/path/to/workflow/
#which workflow (see docs)?
workflow=fltnc
#Which queue runs the workflow job?
workflow_dispatcher_queue=basement
#which nextflow bin?
nextflow_bin=nextflow
#which group?
GROUP="team152"
#####################
#####################

exec_time=$(date -d "today" +"%Y%m%d-%H%M%S")
workdir="${workflow}""${workflow}"
mkdir -p $workdir
#Stub command
MEM=16000; THREADS=2; JOB_NAME="isoseq_dispatcher"; o=dispatcher.o; e=dispatcher.e;
rm $o; rm $e;
cmd="nohup ${nextflow_bin} run ./isoseq2.nf -w ${workdir} -config ./sanger.config -params-file ./localparams.yaml -resume ${workdir} -entry ${workflow} -with-trace -with-report report-${exec_time}.html -with-dag dag-${exec_time}.html > nextflow.nohup.log 2>&1 &"
bsub -R "rusage[mem=${MEM}] select[mem>${MEM}] span[hosts=1]" -q $workflow_dispatcher_queue -M "${MEM}" -n "${THREADS}" -J "${JOB_NAME}" -G "${GROUP}" -o "${o}" -e "${e}" "${cmd}"
echo $cmd

