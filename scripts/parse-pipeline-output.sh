#! /bin/bash

# $1 - directory where annotated_exonic_variants.txt and title_file.txt can be found


# Called Variants
# pull out sample, normal used, chr, ref, alt, mutation, gene, patient id, alt count, total count, ratio for samples in Innovation pools from annotated exonic variants
# add corresponding DMP values in the same row
awk -v titleFile="$1/title_file.txt" 'BEGIN{FS=OFS="\t"; while(getline < titleFile){patient[$3]=$5}}{if(NR==1){for(i=37; i<=NF; i++){sampleIndex[$i]=i} next;} sample=$1; p=patient[sample]; if(index($1, "DummyPool")==0 || seen[p]) next; seen[p]=1; print p, sampleIndex[sample]}' $1/annotated_exonic_variants.txt > PatientID-DMPColumn.txt

# called variants
printf "Sample\tMatchedNormal\tChr\tPosition\tRef\tAlt\tGene\tMutation\tPatient\tAltAlleleReads\tTotalReads\tAltAlleleFraction\tOtherAlleleReads\tDMPRatio\tDMPFraction\n" > called-variants.txt

awk -v titleFile="$1/title_file.txt" 'BEGIN{ FS=OFS="\t"; while(getline < titleFile){patient[$3]=$5} while(getline < "PatientID-DMPColumn.txt"){column[$1]=$2}}{sample=$1; if(index($1, "DummyPool")!=0) next; if(index($2, "CELLFREE")==1) matchedNormal="No"; else matchedNormal="Yes"; if(index(sample, "05500-AJ-P-1")!=0 || patient[sample]==null) next; DMPColumn=column[patient[sample]]; if(DMPColumn==0) {dr="-"; df="-"} else {split($DMPColumn, a, ";"); split(a[1], b, "="); split(a[3], c, "="); split(a[4], d, "="); dr=c[2]"/"b[2]; df=d[2]} print sample, matchedNormal, $3, $4, $5, $6, $8, $14, patient[sample], $27, $25, $28, $25-($26+$27), dr, df}' $1/annotated_exonic_variants.txt >> called-variants.txt


# pull genotypes from various samples from annotated_exonic_variants.txt file
# skip AJ samples, skip mutations with 0 reads, skip standard normals, skip duplicate records
awk -v titleFile="$1/title_file.txt" 'BEGIN{ FS=OFS="\t"; while(getline < titleFile){patient[$3]=$5}}{if(NR==1){for(i=37; i<=NF; i++){sampleName[i]=$i}} else {for(i=37; i<=NF; i++){sample=sampleName[i]; key=$3"\t"$4"\t"$5"\t"$6"\t"sample; if(index(sample, "05500-AJ-P-1")!=0 || patient[sample]==null || seen[key]!=null) continue; seen[key]=key; split($i, a, ";"); split(a[1], b, "="); split(a[2], c, "="); split(a[3], d, "="); if(b[2]==0) continue; print $3, $4, $5, $6, $8, $14, sample, patient[sample], d[2], b[2], d[2]/b[2], b[2]-(c[2]+d[2])}} }' $1/annotated_exonic_variants.txt > pipeline-genotypes-samples.txt

#
