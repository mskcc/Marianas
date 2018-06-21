#!/bin/bash

# process the UMI bams present in the given list. Collapsed fastq files will be created in the current folder.

# $1 - a tab separated file where each line contains path a UMI bam file and path to the corresponding normal pileup file


newProjectDir=`pwd`


mkdir FinalBams
mkdir FinalBams/unfiltered
mkdir FinalBams/simplex-duplex
mkdir FinalBams/duplex


while IFS=$'\t' read -r -a array
do
  bam=${array[0]}
  bam=`readlink -f $bam`
  
  pileup=${array[1]}
  
  # accomodate situations where pileup is not given
  if [[ "$pileup" == "" ]]
  then
    pileup="normal-pileup-not-given"

  else
  	pileup=`readlink -f $pileup`
  fi

  sample=`basename $bam`
  sample=${sample/.bam}

  outputFolder="$sample-folder"
  mkdir $outputFolder

  bsub -eo %J.e -oo %J.o -cwd "$outputFolder" -R "rusage[mem=30]" -R "span[ptile=4]" -n 4 -We 12:00 ~/software/bin/process-umi-bam.sh $bam $pileup

done < $1


#
