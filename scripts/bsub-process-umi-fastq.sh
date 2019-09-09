#! /bin/bash

# process UMI fastqs (put UMI in read name, filter etc.), and create sample folder structure in the current directory. So current directory will be treated as the new project directory.

# $1 - fastq folder that contains all the R1 and R2 fastqs

# $2 - UMI style. Currently supported: loop

# R1 fastqs must have "_R1_" in the name and R2 fastqs must have "_R2_" in the name


# right now assuming there is only one pair of fastqs per sample

for f1 in `ls $1/*_R1_*.fastq.gz`
do
	f2=${f1/_R1_/_R2_}
  	bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=16]" -We 4:00 ~/software/bin/process-umi-fastq.sh $2 $f1 $f2
done











#
