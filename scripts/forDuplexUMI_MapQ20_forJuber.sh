#!/bin/bash
export PATH=/home/pererad1/java8/jdk1.8.0_144/bin:$PATH
input_bam=$1
output_folder=$2
fgbio_jar=/ifs/work/bergerm1/pererad1/ColonWeiserAnalysis/fgbio-0.2.0.jar
scratch_dir=~/scratch_fgbio/

mkdir ${scratch_dir}
mkdir ${output_folder}
#echo "samtools view ${input_bam} |  awk '{print $1 "\t" $3 "\t" $4 "\t" $4+length($10)-1}' > ${output_folder}/readNames.bed"
samtools view ${input_bam} |  awk '{print $1 "\t" $3 "\t" $4 "\t" $4+length($10)-1}' > ${output_folder}/readNames.bed

echo "python /ifs/work/bergerm1/pererad1/ColonWeiserAnalysis/duplexUMI.py ${output_folder}/readNames.bed ${output_folder}/Duplex_UMI_for_readNames.fastq"
python /ifs/work/bergerm1/pererad1/ColonWeiserAnalysis/duplexUMI.py ${output_folder}/readNames.bed ${output_folder}/Duplex_UMI_for_readNames.fastq

echo "java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} AnnotateBamWithUmis -i ${input_bam} -f ${output_folder}/Duplex_UMI_for_readNames.fastq -o ${output_folder}/sample_with_UMI.bam"
java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} AnnotateBamWithUmis -i ${input_bam} -f ${output_folder}/Duplex_UMI_for_readNames.fastq -o ${output_folder}/sample_with_UMI.bam

echo "java -jar ${fgbio_jar} --tmp-dir=${scratch_dir}  SortBam --input=${output_folder}/sample_with_UMI.bam --sort-order=Queryname --output ${output_folder}/sample_with_UMI_sorted.bam"
java  -jar ${fgbio_jar} --tmp-dir=${scratch_dir}  SortBam --input=${output_folder}/sample_with_UMI.bam --sort-order=Queryname --output ${output_folder}/sample_with_UMI_sorted.bam

echo "java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} SetMateInformation -i ${output_folder}/sample_with_UMI_sorted.bam -o ${output_folder}/sample_with_UMI_sorted_mateFixed.bam"
java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} SetMateInformation -i ${output_folder}/sample_with_UMI_sorted.bam -o ${output_folder}/sample_with_UMI_sorted_mateFixed.bam

echo "java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} GroupReadsByUmi -s 'paired' -m 20 -f ${output_folder}/Grouped-mapQ20-histogram -i ${output_folder}/sample_with_UMI_sorted_mateFixed.bam -o ${output_folder}/collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20.bam"
java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} GroupReadsByUmi -s 'paired' -m 20 -f ${output_folder}/Grouped-mapQ20-histogram -i ${output_folder}/sample_with_UMI_sorted_mateFixed.bam -o ${output_folder}/collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20.bam

echo "java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} CallDuplexConsensusReads -i ${output_folder}/collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20.bam -o ${output_folder}/duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20.bam"
java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} CallDuplexConsensusReads -i ${output_folder}/collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20.bam -o ${output_folder}/duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20.bam


echo "java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} FilterConsensusReads -i ${output_folder}/duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20.bam -r /ifs/work/bergerm1/pererad1/impact-GRCh37/Homo_sapiens_assembly19.fasta -M=1 -N=20 -o ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1.bam"
java -jar ${fgbio_jar} --tmp-dir=${scratch_dir} FilterConsensusReads -i ${output_folder}/duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20.bam -r /ifs/work/bergerm1/pererad1/impact-GRCh37/Homo_sapiens_assembly19.fasta -M=1 -N=20 -o ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1.bam

echo "samtools sort -n ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1.bam > ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted.bam"
samtools sort -n ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1.bam > ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted.bam

echo "samtools fastq -1 ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted_R1_.fastq -2 ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted_R2_.fastq ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted.bam"
samtools fastq -1 ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted_R1_.fastq -2 ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted_R2_.fastq ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted.bam

for i in $(ls ${output_folder}/filtered*.fastq); do gzip $i; done

sample_name=$(basename ${input_bam} .bam)
base_dir_name=$(dirname ${output_folder})
mkdir ${base_dir_name}/FastqMapQ20_1-1

mv ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted_R1_.fastq.gz ${base_dir_name}/FastqMapQ20_1-1/FulcrumCollapsed_${sample_name}_mapQ20_1-1_R1_.fastq.gz

mv ${output_folder}/filtered_duplexConsensusReads_collapsed-sample_with_UMI_sorted_mateFixed_paired_mapQ20_1-1_qnameSorted_R2_.fastq.gz ${base_dir_name}/FastqMapQ20_1-1/FulcrumCollapsed_${sample_name}_mapQ20_1-1_R2_.fastq.gz

#rm -rf ${output_folder}
