#!/bin/bash

# Separate simplex-only and unfiltered-only bams. Run it from the collapsing folder in which sample-specific collapsing folders are present

# $1 - a tab separated file where each line contains path a UMI bam file and path to the corresponding normal pileup file (the same file used as input for Marianas collapsing)


newProjectDir=`pwd`


mkdir FinalBams/unfiltered-only
mkdir FinalBams/simplex-only


while IFS=$'\t' read -r -a array
do
  bam=${array[0]}
  bam=`readlink -f $bam`
  sample=`basename $bam`
  sample=${sample/.bam}

  outputFolder="$sample-folder"

  bsub -eo %J.e -oo %J.o -cwd "$outputFolder" -R "rusage[mem=30]" -We 1:00 ~/software/bin/separate-bams.sh $sample

done < $1


#
