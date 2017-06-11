#! /bin/bash

# $1 - R1 fastq file

export PATH=/ifs/work/bergerm1/Innovation/software/metagenomics/kraken-0.10.5-beta:/ifs/work/bergerm1/Innovation/software/metagenomics/KronaTools-2.6.1/bin:/opt/common/CentOS_6/gcc/gcc-4.9.3/bin:/ifs/work/bergerm1/Innovation/software/metagenomics/metaphlan2:/opt/common/CentOS_6/perl/perl-5.22.0/bin:/opt/common/CentOS_6/bowtie2/bowtie2-2.2.4:/opt/common/CentOS_6/python/python-2.7.8/bin:/opt/common/CentOS_6/bedtools/bedtools-2.17.0/:/opt/common/CentOS_6/java/jdk1.8.0_31/bin:/opt/common/CentOS_6/samtools/samtools-1.2:/common/lsf/9.1/linux2.6-glibc2.3-x86_64/etc:/common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:$PATH

export LD_LIBRARY_PATH=/opt/common/CentOS_6/gcc/gcc-4.9.3/lib64/:/common/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/common/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:$LD_LIBRARY_PATH


mpa_dir=/ifs/work/bergerm1/Innovation/software/metagenomics/metaphlan2
minikrakenDB=/ifs/work/bergerm1/Innovation/software/metagenomics/kraken-0.10.5-beta/minikraken-db-20141208

fastq1=$1
sampleName=`basename $fastq1`
sampleName=${sampleName/_IGO_*/}

fastq2=${fastq1/_R1_/_R2_}
outFile=`basename $fastq1`
outFile="${sampleName}_profiled_metagenome.txt"
bowtie2outFile="${sampleName}_bowtie2.txt"

echo -e "$fastq1 $fastq2 $outFile $bowtie2outFile "

# run metaphlan2
#metaphlan2.py $fastq1,$fastq2 --bowtie2out $bowtie2outFile --mpa_pkl ${mpa_dir}/db_v20/mpa_v20_m200.pkl --bowtie2db ${mpa_dir}/db_v20/mpa_v20_m200 --input_type fastq --bt2_ps sensitive-local --min_alignment_len 70 > $outFile

#create Krona Plots
#${mpa_dir}/utils/metaphlan2krona.py -p $outFile -k krona/$sampleName.krona
#ktImportText -o krona/$sampleName.html -n Metagenome krona/$sampleName.krona,$sampleName


# run Kraken
kraken --fastq-input  --gzip-compressed --paired --db $minikrakenDB $fastq1 $fastq2 > $sampleName.kraken
kraken-translate --db $minikrakenDB $sampleName.kraken > $sampleName.kraken.labels
kraken-report --db $minikrakenDB $sampleName.kraken > $sampleName.kraken.report
cat $sampleName.kraken.report | awk '$3>=5' > kraken-short/$sampleName.short.txt
echo -e "Innovation Lab Kraken Run" > $sampleName.t
echo -e "Sample: $sampleName" >> $sampleName.t
printf "PercentReads\t#CladeReads\t#TaxonReads\tRank\tNCBI-ID\tName\n" >> $sampleName.t
cat $sampleName.t kraken-short/$sampleName.short.txt > $sampleName.t1
mv $sampleName.t1 kraken-short/$sampleName.short.txt
rm $sampleName.t
kraken-mpa-report --db $minikrakenDB $sampleName.kraken > $sampleName.kraken.mpa

#create Krona Plots
${mpa_dir}/utils/metaphlan2krona.py -p $sampleName.kraken.mpa -k kraken-krona/$sampleName.krona
ktImportText -o kraken-krona/$sampleName.html -n Metagenome kraken-krona/$sampleName.krona,$sampleName




#
