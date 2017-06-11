#! /bin/bash

# run it from the directory where output is desired

# $1 - directory containing bam files

mkdir metaphlan2-krona
mkdir kraken-krona
mkdir kraken-short

for bam in `ls $1/*.bam`
do
  bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=16]" -We 0:60 $HOME/software/bin/metagenomics-from-bams.sh $bam

done















#
