#!/bin/bash

in_f=$1
bam_folder=$2
base_dir=$3
mkdir ${base_dir}
for line in $(cat ${in_f}); do
    base=$(basename ${line} .bam)
    output_f=${base_dir}/${base}
    echo ${base}_Submitting
    bsub -o ${base_dir}/${base}_o.log -e ${base_dir}/${base}_e.log -R "rusage[mem=32]" ~/software/bin/forDuplexUMI_MapQ20_forJuber.sh  ${bam_folder}/${line} ${output_f}
done
