#! /bin/bash

# get the read count for each barcode present in the given fastq file. Also make a separate file that lists read counts for standard barcodes only.

# $1 - fastq.gz file to generate counts for
# $2 - name stub to use


echo -e "counting barcodes..."
zcat $1 | awk '{if(NR%10000000==0){print NR > "/dev/stderr"} if(NR%4!=1) next; split($0, a, ":"); barcode=a[length(a)]; barcodes[barcode]++}END{for(barcode in barcodes) print barcode"\t"barcodes[barcode]}' | sort -k 2,2nr > $2-barcodes.txt

echo -e "separating standard barcodes..."
awk 'BEGIN{FS=OFS="\t"; while(getline < "barcodeKey96.txt"){barcodes[$1]=$2}}{if(barcodes[$1]!=null){print $1, $2, barcodes[$1]}}' $2-barcodes.txt > $2-standard-barcodes.txt

#
