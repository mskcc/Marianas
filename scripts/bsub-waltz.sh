#! /bin/bash

# $1 - directory containing bam files
# $2 - module to run: Metrics, Genotyping, SignatureFinding
# $3 - minimum mapping quality
# $4 - reference fasta file
# $5 - intervals bed file
# $6 - module argument
	# empty for Metrics module
	# loci bed file for Genotyping module
	# comma-separated list of signatures to look for for SignatureFinding module

folder=$1

for f in `ls $folder/*.bam`
do
	#bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=8]" -We 0:60 ~/software/bin/locuspocus.sh $2 $3 $f $4 $5 $6
	bsub -eo %J.e -oo %J.o -cwd "." -R "rusage[mem=8]" -We 0:60 ~/software/bin/waltz.sh $2 $3 $f $4 $5 $6
done
