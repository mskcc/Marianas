#! /bin/bash

printf "Sample\tChr\tPosition\tRef\tTotal\tA\tC\tG\tT\n" > FP-allele-fractions.txt

# collect numbers from Waltz pileup
#for f in `ls *-pileup-without-duplicates.txt`
for f in `ls *-pileup.txt`
do
  sample="${f/_bc*/}"
  #sample=`echo $sample | awk 'BEGIN{FS="_bc"}{print $1}'`

  # add FP allele fractions to the file
  awk -v sample=$sample 'BEGIN{FS=OFS="\t"; while(getline < "/Users/patelj1/workspace/Waltz/MutationsFiles/Fingerprinting-SNPs.txt"){key=$1"\t"$2; pos[key]=key}}{key=$1"\t"$2; if(pos[key]!=null) {total=$5+$6+$7+$8; if(total==0) {total=1}; print sample, $1, $2, $3, $4, $5/total, $6/total, $7/total, $8/total}}' $f >> FP-allele-fractions.txt

done




#
