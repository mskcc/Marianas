#! /bin/bash

# process UMI fastqs (put UMI in read name, filter etc.), and create sample folder structure in the current directory. So current directory will be treated as the new project directory.

# $1 - project folder that contains sample folders that in turn contain fastq files and sample sheets

# $2 - UMI style. Currently supported: loeb and loop

newProjectDir=`pwd`

# right now assuming there is only one pair of fastqs per sample

for f1 in `ls $1/*/*_R1_*.fastq.gz`
do
  bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=16]" -We 4:00 ~/software/bin/process-umi-fastq.sh $2 $newProjectDir $f1
done











#
