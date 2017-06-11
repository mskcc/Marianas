#! /bin/bash

# generate the bed file that is valid to use for Waltz genotyping for the given mutations file. Still might require manual inspection and modification.

# $1 - mutations file (must end in .txt)




awk '{FS=OFS="\t"}{width=length($3)+length($4)+10; start=$2-width; end=$2+width; print $1, start, end, "-", $5}' $1 | sort -k 2,2n > t.bed

# merge consecutive intervals that are less than 30bp apart
awk 'BEGIN{FS=OFS="\t"; dist=30}{if($1!=chr){if(chr!=null){print chr, start, stop, "-", name;} chr=$1; start=$2; stop=$3; name=$5; next} if($2-stop<dist && $2-stop>=0){ stop=$3;} else if(start-$3<dist && start-$3>=0){start=$2} else{print chr, start, stop, "-", name; chr=$1; start=$2; stop=$3; name=$5;} }END{print chr, start, stop, "-", name;}' t.bed > ${1/.txt/.bed}

rm t.bed



#
