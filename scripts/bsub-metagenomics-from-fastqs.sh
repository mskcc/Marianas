#! /bin/bash

# run it from the directory where output is desired

# $1 - Project fastq directory

mkdir metaphlan2-krona
mkdir kraken-krona
mkdir kraken-short

for fastq1 in `ls $1/Sample_*/*_R1_*.fastq.gz`
do
  bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=16]" -We 0:60 $HOME/software/bin/metagenomics-from-fastqs.sh $fastq1

done















#
