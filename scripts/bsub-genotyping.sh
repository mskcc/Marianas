#! /bin/bash

# $1 - directory containing bam files
# $2 - intervals bed file
# $3 - mutations file

referenceFasta=~/resources/impact-GRCh37/Homo_sapiens_assembly19.fasta

folder=$1

for f in `ls $folder/*.bam`
do
	bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=8]" -We 0:60 ~/software/bin/waltz.sh Genotyping 20 $f $referenceFasta $2 $3
done
