#!/bin/bash

# $1 - Collapsing run directory i.e. parent directory of collapsed bam directories. Must have on-target-positions.txt in the current directory


# process coverage
printf "Sample\tDuplexFamilies\tDuplexFamiliesWithoutSingletons\n" > duplex-families.txt


# collect mean coverage from Waltz Metrics output
for f in `ls -d $1/*-folder`
do
	sample=`basename $f`
	sample="${sample/-IGO*/}"
	printf "$sample\t" >> duplex-families.txt

	numbers=`awk 'BEGIN{FS=OFS="\t"; while(getline < "on-target-positions.txt"){positions[$1"\t"$2]=$1"\t"$2}}{if(positions[$1"\t"$2]==null) next; total++; if($4+$5>1){total2++} if($4!=0 && $5 !=0){ds++}}END{print ds/total, ds/total2}' $f/first-pass.txt`

	#numbers=`awk 'BEGIN{FS=OFS="\t"; while(getline < "on-target-positions.txt"){positions[$1"\t"$2]=$1"\t"$2}}{if(positions[$1"\t"$2]==null) next; total++; if($4!=0 && $5 !=0){ds++} if($4>1 || $5>1){total2++; if($4!=0 && $5 !=0){ds2++}}}END{print ds/total, ds2/total2}' $f/first-pass.txt` 

	printf "$numbers\n" >> duplex-families.txt

done
