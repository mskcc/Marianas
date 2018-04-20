#! /bin/bash

# get the UMI cluster details at the given position for the given bam

# $1 - position (eg. 3:14455763-14455763, samtools-view style)

# $2 - original (uncollapsed) bam file

# $3 - pileup file corresponsing to the bam file (made with correct bed file)

export TMPDIR=/ifs/work/scratch

java=/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java
referenceFasta=~/resources/impact-GRCh37/Homo_sapiens_assembly19.fasta
bamFile=`basename $2`
bamFile=${bamFile/.bam/}-${1/:/-}.bam
detailsFile=${bamFile/.bam/}-cluster-details.txt

echo -e "making bam of reads that cover $1"
samtools view -bh $2 $1 > $bamFile

echo -e "gathering cluster details"
$java -server -Xms8g -Xmx8g -cp ~/software/Marianas.jar org.mskcc.marianas.umi.duplex.debug.ClusterDetailsAtPosition $bamFile $3 $1 1 2 $referenceFasta $detailsFile
