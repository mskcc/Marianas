#!/bin/bash

# run MuTect

# $1 - tumor bam
# $2 - normal bam
# $3 - mutect run name (used for output files)

#java=/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java
java=/opt/common/CentOS_6/java/jdk1.6.0_45/bin/java
referenceFasta=~/resources/impact-GRCh37/Homo_sapiens_assembly19.fasta
mutectJar=/opt/common/CentOS_6/muTect/muTect-1.1.4/muTect.jar
dbsnp=/ifs/e63data/bergerm1/Resources/DMP/pubdata/dbSNP/VERSIONS/dbsnp_v137/dbsnp_137.b37.vcf
cosmic=/ifs/e63data/bergerm1/Resources/DMP/pubdata/cosmic/VERSIONS/cosmic_v54_120711/hg19_cosmic_v54_120711.vcf

export TMPDIR=/ifs/work/scratch

# run MuTect
 $java -server -Xms4g -Xmx8g -jar $mutectJar -T MuTect --input_file:normal $2 --input_file:tumor $1 --reference_sequence $referenceFasta --dbsnp $dbsnp --cosmic $cosmic -o $3-mutect.txt -vcf $3-mutect.vcf --enable_extended_output  -dcov 50000 -rf BadCigar -rf MappingQuality -mmq 20



#
