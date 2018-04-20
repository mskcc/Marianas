#!/bin/bash

# call variants using Marianas variant caller

# $1 - a tab separated file where each line contains path to tumor pileup file, normal pileup file and a sample name

hotspotsFile=/home/patelju1/resources/hotspots/hotspot-list-union-v1-v2.txt
noiseFrequenciesFile=/home/patelju1/resources/polishing-normals-marianas-3-3-v2/noise-frequencies.txt


newProjectDir=`pwd`

while IFS=$'\t' read -r -a array
do
  tumorPileup=${array[0]}
  normalPileup=${array[1]}
  sample=${array[2]}
  bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=16]" -We 0:59 ~/software/bin/marianas-variant-calling.sh `readlink -f $tumorPileup` `readlink -f $normalPileup` $sample $hotspotsFile $noiseFrequenciesFile

done < $1


#
