#! /bin/bash

# $1 - directory containing bam files
# $2 - bed file of targeted regions

folder=$1
bedFile=$2

coverageThreshold="50"
geneList="/home/patelju1/resources/gene-list/juber-hg19-gene-list.bed"

for f in `ls $folder/*.bam`
do
	bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=8]" -We 6:00 ~/software/bin/bam-metrics.sh $f $coverageThreshold $geneList $bedFile
done
