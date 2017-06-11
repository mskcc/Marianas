#! /bin/bash

# $1 - title file, optional. It is used to get patient ID.


printf "Sample\tPatient\tChr\tPosition\tRef\tAlt\tEventType\tEventName\tTotalSupportingReads\tTotalReads\tTotalAltAlleleFreq\tUniqueSupportingReads\tUniqueReads\tUniqueAltAlleleFreq\n" > genotypes.txt

# collect numbers from Waltz genotyping
for f in `ls *-genotypes.txt`
do
  sample="${f/-genotypes.txt/}"
  sample=`echo $sample | awk 'BEGIN{FS="_bc"}{print $1}'`

  # add genotypes to genotypes.txt
	awk -v sample=$sample -v titleFile=$1 'BEGIN{FS=OFS="\t"; if(titleFile!=null){while(getline < titleFile){patients[$3]=$5}}}{if($8==0) tFraction=0; else tFraction=$7/$8; if($10==0) uFraction=0; else uFraction=$9/$10; patient=patients[sample]; if(patient==null) patient="-"; print sample, patient, $1, $2, $3, $4, $5, $6, $7, $8, tFraction, $9, $10, uFraction}' $f >> genotypes.txt

done




#
