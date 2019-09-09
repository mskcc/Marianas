#! /bin/bash

# align paired end fastq files, sort sam, create bams and indexes, mark duplicates, fix mate information. Run it from where the bam files are desired.

# $1 - directory where fastq files are located.

# $2 - run PICARD MarkDuplicates. yes or no. Default is yes.


for f in `ls $1/*_R1_*.fastq.gz`
do
  bsub -q sol -cwd "." -R "rusage[mem=30]" -R "span[ptile=4]" -n 4 -We 1:00 -o %J.o -e %J.e ~/software/bin/juber-fastq-to-bam.sh $f $2

done



#
