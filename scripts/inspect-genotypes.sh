#! /bin/bash

# $1 - directory containing bam files
# $2 - intervals bed file
# $3 - mutations file

referenceFasta=~/resources/impact-GRCh37/Homo_sapiens_assembly19.fasta

tempDir=`mktemp -d -t inspect-genotypes.XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX`
cd $tempDir

folder=$1

printf "Sample\tChr\tPosition\tRef\tAlt\tEventType\tEventName\tTotalSupportingReads\tTotalReads\tTotalAltAlleleFreq\tUniqueSupportingReads\tUniqueReads\tUniqueAltAlleleFreq\n" > genotypes.txt

for f in `ls $folder/*.bam`
do
	~/software/bin/waltz.sh Genotyping QualityFilter $f $referenceFasta $2 $3

	baseName=`basename $f`
	genotypesFile=${baseName/.bam/-genotypes.txt}
	sample="${genotypesFile/-genotypes.txt/}"
  sample=`echo $sample | awk 'BEGIN{FS="_bc"}{print $1}'`

	echo -e $f

	if [ ! -e  "$genotypesFile" ]
	then
		continue
	fi

	cat $genotypesFile

	# add genotypes to genotypes.txt
	awk -v sample=$sample 'BEGIN{FS=OFS="\t"}{if($8==0) tFraction=0; else tFraction=$7/$8; if($10==0) uFraction=0; else uFraction=$9/$10; print sample, $1, $2, $3, $4, $5, $6, $7, $8, tFraction, $9, $10, uFraction}' $genotypesFile >> genotypes.txt
done

echo -e "temp dir: $tempDir"






#
