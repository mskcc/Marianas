#! /bin/bash

# $1 - bam file aligned to human genome

export PATH=/ifs/work/bergerm1/Innovation/software/metagenomics/kraken-0.10.5-beta:/ifs/work/bergerm1/Innovation/software/metagenomics/KronaTools-2.6.1/bin:/opt/common/CentOS_6/gcc/gcc-4.9.3/bin:/ifs/work/bergerm1/Innovation/software/metagenomics/metaphlan2:/opt/common/CentOS_6/perl/perl-5.22.0/bin:/opt/common/CentOS_6/bowtie2/bowtie2-2.2.4:/opt/common/CentOS_6/python/python-2.7.8/bin:/opt/common/CentOS_6/bedtools/bedtools-2.17.0/:/opt/common/CentOS_6/java/jdk1.8.0_31/bin:/opt/common/CentOS_6/samtools/samtools-1.2:/common/lsf/9.1/linux2.6-glibc2.3-x86_64/etc:/common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:$PATH

export LD_LIBRARY_PATH=/opt/common/CentOS_6/gcc/gcc-4.9.3/lib64/:/common/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/common/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:$LD_LIBRARY_PATH


mpa_dir=/ifs/work/bergerm1/Innovation/software/metagenomics/metaphlan2
minikrakenDB=/ifs/work/bergerm1/Innovation/software/metagenomics/kraken-0.10.5-beta/minikraken-db-20141208


bam=$1
sampleName=`basename $bam`
sampleName=${sampleName/-IGO-*/}
fastq=${sampleName}.fastq
outFile="${sampleName}_profiled_metagenome.txt"
bowtie2outFile="${sampleName}_bowtie2.txt"

echo -e "$bam $fastq $outFile $bowtie2outFile "

#extract fastq of unmapped reads
samtools view -f 4 $bam | awk 'BEGIN{FS=OFS="\t"}{print "@"$1; print $10; print "+"; print $11}' > $fastq

# run MetaPhlAn2
# metaphlan2.py $fastq --bowtie2out $bowtie2outFile --mpa_pkl ${mpa_dir}/db_v20/mpa_v20_m200.pkl --bowtie2db ${mpa_dir}/db_v20/mpa_v20_m200 --input_type fastq --bt2_ps very-sensitive-local --min_alignment_len 30 > $outFile
metaphlan2.py $fastq --bowtie2out $bowtie2outFile --mpa_pkl ${mpa_dir}/db_v20/mpa_v20_m200.pkl --bowtie2db ${mpa_dir}/db_v20/mpa_v20_m200 --input_type fastq --bt2_ps sensitive-local --min_alignment_len 70 > $outFile

#create Krona Plots
${mpa_dir}/utils/metaphlan2krona.py -p $outFile -k metaphlan2-krona/$sampleName.krona
ktImportText -o metaphlan2-krona/$sampleName.html -n Metagenome metaphlan2-krona/$sampleName.krona,$sampleName


# run Kraken
kraken --fastq-input  --db $minikrakenDB $fastq > $sampleName.kraken
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
