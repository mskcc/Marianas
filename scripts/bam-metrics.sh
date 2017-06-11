#! /bin/bash

# $1 - bam file
# $2 - coverage threshold to apply when computing covered regions
# $3 - gene list to use to annoatate covered regions
# $4 - bed file of targeted regions

java=/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java
referenceFasta=~/resources/impact-GRCh37/Homo_sapiens_assembly19.fasta

bamFile=$1
coverageThreshold=$2
geneList=$3
bedFile=$4

export TMPDIR=/ifs/work/scratch

# scan the bam
$java -server -Xms4g -Xmx4g -cp ~/software/Marianas.jar org.mskcc.marianas.metrics.CountReads $bamFile $coverageThreshold $geneList $bedFile

# get coverage and other metrics from Waltz
$java -server -Xms4g -Xmx4g -cp ~/software/Waltz.jar org.mskcc.juber.waltz.Waltz Metrics 20 $bamFile $referenceFasta $bedFile
