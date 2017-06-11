#! /bin/bash


# create copy number files
printf "Interval\tSample\tTotalCoverage\tGC\n" > waltz-cn-total.txt
printf "Interval\tSample\tUniqueCoverage\tGC\n" > waltz-cn-unique.txt


for f in `ls *-intervals.txt`
do
  sample="${f/-intervals.txt/}"
  sample=`echo $sample | awk 'BEGIN{FS="_bc"}{print $1}'`

  # do intra-sample normalization
  # use average cn probe coverage as normalizing factor
  normFactor=`awk 'BEGIN{FS=OFS="\t"; count=0;}{if(index($4, "cn_")!=1) next; count++; coverage += $7}END{print coverage/count}' $f`

  awk -v sample=$sample -v normFactor=$normFactor 'BEGIN{FS=OFS="\t"}{if(index($4, "FP_")==1) next; split($4, a, "_"); print $4, sample, $7/normFactor, $8;}' $f >> waltz-cn-total.txt


  # sample mean as the normalizing factor
  normFactor=`awk 'BEGIN{FS=OFS="\t"; count=0;}{if(index($4, "cn_")!=1) next; count++; coverage += $7}END{print coverage/count}' ${f/-intervals/-intervals-without-duplicates}`

  awk -v sample=$sample -v normFactor=$normFactor 'BEGIN{FS=OFS="\t";}{if(index($4, "FP_")==1) next; split($4, a, "_"); print $4, sample, $7/normFactor, $8;}'  ${f/-intervals/-intervals-without-duplicates} >> waltz-cn-unique.txt
done










####
