#! /bin/bash

# make paired end fastq files, one for each given barcode, from given pair of fastq files

# $1 - R1 fastq.gz file from which to make smaller fastq files
# $2 - tab-delimited file of barcode, read count for that barcode and barcode name. R1 and R2 files will be generated for each barcode.

id=`basename $2`
id=${id/.txt/}

while IFS=$'\t' read -r -a array
do
  barcode=${array[0]}
  count=${array[1]}
  barcodeName=${array[2]}
  R1Name="$id-$barcodeName-$barcode-$count-R1.fastq"
  R2Name="$id-$barcodeName-$barcode-$count-R2.fastq"

  echo -e "making fastq pair for $barcode"
  zcat $1 | awk -v barcode="$barcode" '{getline seq; getline garbage; getline qual; split($0, a, ":"); bc=a[length(a)]; if(bc==barcode) {print $0; print seq; print garbage; print qual}}' > $R1Name

  zcat ${1/_R1_/_R2_} | awk -v barcode="$barcode" '{getline seq; getline garbage; getline qual; split($0, a, ":"); bc=a[length(a)]; if(bc==barcode) {print $0; print seq; print garbage; print qual}}' > $R2Name

  # zip 'em
  gzip $R1Name
  gzip $R2Name

done < $2




#
