#!/bin/bash

# Process the given UMI bam file and generate collapsed fatsqs and stats

# $1 - UMI bam file to be processed
# $2 - normal pileup file corresponding to the bam file (made with correct bed file)


export TMPDIR=/ifs/work/scratch


java="/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java"
referenceFasta=~/resources/impact-GRCh37/Homo_sapiens_assembly19.fasta
#referenceFasta=/ifs/work/bergerm1/Innovation/projects/Juber/Wendy-Viral/resources/hg19-mcpyv-ebv-hpv/hg19-mcpyv-ebv-hpv.fasta

minMappingQuality="1"
minBaseQuality="20"
mismatches="0"
wobble="1"
minConsensusPercent="90"
jarFile="$HOME/software/Marianas.jar"


echo -e "minMappingQuality=$minMappingQuality, minBaseQuality=$minBaseQuality, mismatches=$mismatches, wobble=$wobble, minConsensusPercent=$minConsensusPercent"

bam=$1
pileup=$2
sample=`basename $bam`
sample=${sample/.bam}


# pass 1
echo -e "`date` Pass 1"
$java -server -Xms8g -Xmx8g -cp $jarFile org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqFirstPass $bam $pileup $minMappingQuality $minBaseQuality $mismatches $wobble $minConsensusPercent $referenceFasta


# sort the first pass output by mate position
echo -e "`date` Sorting first-pass.txt"
sort -S 8G -k 6,6n -k 8,8n first-pass.txt > first-pass.mate-position-sorted.txt


# pass 2
echo -e "`date` Pass 2"
$java -server -Xms8g -Xmx8g -cp $jarFile org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqSecondPass $bam $pileup $minMappingQuality $minBaseQuality $mismatches $wobble $minConsensusPercent $referenceFasta first-pass.mate-position-sorted.txt

# delete unnecessory files
# rm first-pass.mate-position-sorted.txt

# gzip
echo -e "`date` Compressing fastqs"
gzip collapsed_R1_.fastq
gzip collapsed_R2_.fastq


# make collapsed bam
echo -e "`date` Running juber-fastq-to-bam.sh "
~/software/bin/juber-fastq-to-bam.sh collapsed_R1_.fastq.gz no


# separate duplex and simplex-duplex bams
echo -e "`date` Separating simplex-duplex and duplex bams "
$java -server -Xms8g -Xmx8g -cp $jarFile org.mskcc.marianas.umi.duplex.postprocessing.SeparateBams collapsed.bam



# link in the FinalBams folder
echo -e "`date` Linking bams"
wd=`readlink -f .`

cd ../FinalBams/unfiltered
ln -s $wd/collapsed.bam $sample-unfiltered.bam
ln -s $wd/collapsed.bai $sample-unfiltered.bai
ln -s $wd/collapsed.bai $sample-unfiltered.bam.bai

cd ../simplex
ln -s $wd/collapsed-simplex.bam $sample-simplex.bam
ln -s $wd/collapsed-simplex.bai $sample-simplex.bai
ln -s $wd/collapsed-simplex.bai $sample-simplex.bam.bai

cd ../duplex
ln -s $wd/collapsed-duplex.bam $sample-duplex.bam
ln -s $wd/collapsed-duplex.bai $sample-duplex.bai
ln -s $wd/collapsed-duplex.bai $sample-duplex.bam.bai




echo -e "`date` Done."





#
