#! /bin/bash

# process the UMI bams present in the given folder. Collapsed fastq files will be created in the current folder

# $1 - a tab separated file where each line contains path a UMI bam file and path to the corresponding normal pileup file

newProjectDir=`pwd`


mkdir FinalBams

while IFS=$'\t' read -r -a array
do
  bam=${array[0]}
  pileup=${array[1]}
  sample=`basename $bam`
  sample=${sample/.bam}
  outputFolder="$sample-folder"
  mkdir $outputFolder
  bsub -eo %J.e -oo %J.o -cwd "$outputFolder" -R "rusage[mem=30]" -n 4 -We 12:00 ~/software/bin/process-umi-bam.sh `readlink -f $bam` `readlink -f $pileup`

done < $1


#
