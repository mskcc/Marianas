#! /bin/bash

# run in MantaOutput directory

# $1 - pairing file with tumor-normal pairs


printf "Tumor\tNormal\tChr\tPosition\tRef\tAlt\tSVTYPE\tBND_DEPTH\tMATE_BND_DEPTH\tSOMATICSCORE\tJUNCTION_SOMATICSCORE\trefPairedReads\trefSplitReads\taltPairedReads\taltSplitReads\n" > manta-output.txt

# collect output from samples
while IFS=$'\t' read -r -a array
do
  tumor=${array[0]}
  normal=${array[1]}

  #gunzip $tumor/result/variants/somaticSV.vcf.gz
  vcf=$tumor/results/variants/somaticSV.vcf

  # add SVs to the aggregate file
  awk -v tumor="$tumor" -v normal="$normal" 'BEGIN{FS=OFS="\t"}{if(index($1, "#")==1) next; delete a; delete m; n=split($8, a, ";"); for(i=1; i<=n; i++){delete b; split(a[i], b, "="); m[b[1]]=b[2]; } delete a; split($11, a, ":"); split(a[1], b, ","); refP=b[1]; altP=b[2]; split(a[2], b, ","); refS=b[1]; altS=b[2];   print tumor, normal, $1, $2, $4, $5, m["SVTYPE"], m["BND_DEPTH"], m["MATE_BND_DEPTH"], m["SOMATICSCORE"], m["JUNCTION_SOMATICSCORE"], refP, refS, altP, altS}' $vcf >> manta-output.txt

done < $1




#
