#!/bin/bash

# Separate simplex-only and unfiltered-only bams

# $1 - sample name

sample=$1


export TMPDIR=/ifs/work/scratch


java="/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java"
jarFile="$HOME/software/Marianas-1.5.0-sb.jar"



# separate duplex and simplex bams
echo -e "`date` Separating simplex-only and unfiltered-only bams "
$java -server -Xms8g -Xmx8g -cp $jarFile org.mskcc.marianas.umi.duplex.postprocessing.SeparateBamsSimplexOnly collapsed.bam


#touch collapsed-simplex-only.bam
#touch collapsed-unfiltered-only.bam

#touch collapsed-simplex-only.bai
#touch collapsed-unfiltered-only.bai


# link in the FinalBams folder
echo -e "`date` Linking bams"
wd=`readlink -f .`

cd ../FinalBams/unfiltered-only
ln -s $wd/collapsed-unfiltered-only.bam $sample-unfiltered-only.bam
ln -s $wd/collapsed-unfiltered-only.bai $sample-unfiltered-only.bai
ln -s $wd/collapsed-unfiltered-only.bai $sample-unfiltered-only.bam.bai

cd ../simplex-only
ln -s $wd/collapsed-simplex-only.bam $sample-simplex-only.bam
ln -s $wd/collapsed-simplex-only.bai $sample-simplex-only.bai
ln -s $wd/collapsed-simplex-only.bai $sample-simplex-only.bam.bai


echo -e "`date` Done."





#
