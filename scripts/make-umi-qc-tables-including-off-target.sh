#!/bin/bash

# $1 - The collapsed bam folder that contains sample folders that contain collapsed bams, collapsed fastqs and first-pass.txt. 


# make output files
echo -e "ClusterSize\tSample" > cluster-sizes.txt
echo -e "ClusterSize\tSample" > cluster-sizes-post-filtering.txt
echo -e "Clusters\tSample" > clusters-per-position.txt
echo -e "Clusters\tSample" > clusters-per-position-post-filtering.txt
echo -e "Sample\tDuplexFamilies\tDuplexFamiliesWithoutSingletons" > duplex-families.txt


# process samples
for sampleFolder in `ls -d $1/*-folder`
do
  sampleName=`basename $sampleFolder`
  sampleName=${sampleName/_IGO*/}

  # cluster sizes
  awk -v sample=$sampleName 'BEGIN{FS=OFS="\t";}{print $4+$5, sample}' $sampleFolder/first-pass.txt >> cluster-sizes.txt

  # cluster sizes post filtering
  gunzip -c $sampleFolder/collapsed_R1_.fastq.gz | grep "@Marianas" | awk -v sample=$sampleName 'BEGIN{FS=OFS="\t"; FS=":"}{print $5+$6, sample}' >> cluster-sizes-post-filtering.txt

  # clusters per position
  awk -v sample=$sampleName 'BEGIN{FS=OFS="\t";}{print $1"AND"$2}' $sampleFolder/first-pass.txt | sort | uniq -c | awk -v sample=$sampleName '{print $1"\t"sample}' >> clusters-per-position.txt

  # clusters per position post filtering
  gunzip -c $sampleFolder/collapsed_R1_.fastq.gz | grep "@Marianas" | awk -v sample=$sampleName 'BEGIN{FS=OFS="\t"; FS=":"}{print $3"AND"$4}' | sort | uniq -c | awk -v sample=$sampleName '{print $1"\t"sample}' >> clusters-per-position-post-filtering.txt

  # duplex families
  printf "$sampleName\t" >> duplex-families.txt
  numbers=`awk 'BEGIN{FS=OFS="\t";}{total++; if($4+$5>1){total2++} if($4!=0 && $5 !=0){ds++}}END{print ds/total, ds/total2}' $sampleFolder/first-pass.txt`
  printf "$numbers\n" >> duplex-families.txt

done




#
