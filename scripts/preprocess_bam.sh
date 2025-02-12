#!/usr/bin/env bash

sample_id=$1
bam=$2
chrom=$3
suffix=$4
cpus=$5




output_bam="${sample_id}.${chrom}${suffix}.bam"
#Removes supplementary alignments (note this may affect the detection of chimeric/fusion genes)
#Appends sample_id to CB tag
#Adds sample_id to qname (read id)

samtools view -@ ${cpus} -h -F 2048 ${bam} ${chrom} |\
awk -v sample=${sample_id} -v chrom=${chrom} 'BEGIN {OFS="\t"}{\
if ($0 ~ /^@/) {\
    if($0 ~ /^@SQ/) {if ($2=="SN:"chrom) {print} else{next}} else{print};\
next} \
for(i=12; i<=NF; i++) \
  {if ($i ~ /^CB:Z:/) {$i = $i"_"sample}}print}' |\
awk -v sample=${sample_id} 'BEGIN {OFS="\t"}{if ($0 ~ /^@/) {print; next} else {$1 = $1"_"sample;print $0}}' |\
samtools view -@ ${cpus} -h -bo $output_bam - ;
samtools index -@ ${cpus} $output_bam ;
