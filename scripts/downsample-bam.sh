#!/bin/bash

# $1 - bam file
# $2 - fraction to keep

seed=54661
out=`basename $1`

echo -e  "-s is: $seed.${2/*./}"
echo -e "downsampling..."
samtools view -bh -s $seed.${2/*./} $1 > $out

echo -e "sorting..."
samtools sort $out -o $out.sorted
mv $out.sorted $out

# mark duplicates again! This reomves incorrect duplicate flags and adds new duplicate flags
echo -e "marking duplicates..."
/opt/common/CentOS_6/java/jdk1.7.0_75/bin/java -Xmx16g -jar /opt/common/CentOS_6/picard/picard-tools-1.96//MarkDuplicates.jar I=$out O=${out/.bam/-MD.bam} ASSUME_SORTED=true METRICS_FILE=${out}-MD.metrics TMP_DIR=/ifs/work/scratch/ COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

mv ${out/.bam/-MD.bam} $out
mv ${out/.bam/-MD.bai} ${out/.bam/.bai}

echo -e "done."

