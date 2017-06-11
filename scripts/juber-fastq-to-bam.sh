#! /bin/bash

# $1 - R1 fastq. Must have _R1_ in the name and must end in .fastq.gz. Corresspodning R2 file must exist in the same directory

# $2 - run PICARD MarkDuplicates. yes or no. Default is yes.

export TMPDIR=/ifs/work/scratch

workingDir=`pwd`
# referenceFasta=/ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta
referenceFasta=/ifs/work/bergerm1/Innovation/projects/Juber/Wendy-Viral/resources/hg19-mcpyv-ebv-hpv/hg19-mcpyv-ebv-hpv.fasta
R1Fastq=$1
R2Fastq=${R1Fastq/_R1_/_R2_}
R1FastqBasename=`basename $R1Fastq`
R2FastqBasename=`basename $R2Fastq`
sample=${R1FastqBasename/_R1_*/}
barcode=`gunzip -c $R1Fastq | head -1 | awk '{split($0, a, ":"); print a[length(a)]}'`

# clip barcode residues
echo -e "`date` clipping adapter residues"
/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl /opt/common/CentOS_6/trim_galore/Trim_Galore_v0.2.5/trim_galore --paired --gzip -q 1 --suppress_warn --stringency 3 -length 25 -o $workingDir -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGATTCATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT $R1Fastq $R2Fastq

# align fastqs
echo -e "`date` aligning fastqs"
/opt/common/CentOS_6/bwa/bwa-0.7.5a/bwa mem -t 4 -PM -R "@RG\tID:${sample}\tLB:Garbage\tSM:${sample}\tPL:Illumina\tPU:${barcode}\tCN:InnovationLab" $referenceFasta ${R1FastqBasename/.fastq.gz/_cl.fastq.gz} ${R2FastqBasename/.fastq.gz/_cl.fastq.gz} > ${sample}.sam

rm ${R1FastqBasename/.fastq.gz/_cl.fastq.gz} ${R2FastqBasename/.fastq.gz/_cl.fastq.gz}

# sort sam and make bam
echo -e "`date` sorting sam and making bam"
/opt/common/CentOS_6/java/jdk1.7.0_75/bin/java -Xmx16g -jar /opt/common/CentOS_6/picard/picard-tools-1.96//AddOrReplaceReadGroups.jar I=${sample}.sam O=${sample}.bam SO=coordinate RGID=${sample} RGLB=Garbage RGPL=Illumina RGPU=${barcode} RGSM=${sample} RGCN=InnovationLab TMP_DIR=/ifs/work/scratch/ COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

rm ${sample}.sam

if [ "$2" == "no" ]
then
  echo -e "`date` NOT marking duplicates!"
else
  # mark duplicates
  echo -e "`date` marking duplicates"
  /opt/common/CentOS_6/java/jdk1.7.0_75/bin/java -Xmx16g -jar /opt/common/CentOS_6/picard/picard-tools-1.96//MarkDuplicates.jar I=${sample}.bam O=${sample}_aln_srt_MD.bam ASSUME_SORTED=true METRICS_FILE=${sample}_aln_srt_MD.metrics TMP_DIR=/ifs/work/scratch/ COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

  mv ${sample}_aln_srt_MD.bam ${sample}.bam
  mv ${sample}_aln_srt_MD.bai ${sample}.bai
fi


# fix mate information
echo -e "`date` fixing mate information"
/opt/common/CentOS_6/java/jdk1.7.0_75/bin/java -Xmx24g -jar /opt/common/CentOS_6/picard/picard-tools-1.96//FixMateInformation.jar I=${sample}.bam O=${sample}_FX.bam SO=coordinate TMP_DIR=/ifs/work/scratch/ COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT


mv ${sample}_FX.bam ${sample}.bam
mv ${sample}_FX.bai ${sample}.bai
ln -s ${sample}.bai ${sample}.bam.bai




#
