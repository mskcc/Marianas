#! /bin/bash

# $1 - file with gene names on a single line, separated by space. Do not use tabs. Just the gene name (eg. ROS1) picks up all its exons. Gene name followed by "_introns" (eg. ROS1_introns) picks up all its introns. All fingerprinting and copy number probes are included automatically. Unnecessary probes should be removed manually.

masterListFile=~/resources/IDT-probes/batch-2.txt

genes="`head -1 $1`"

# add introns and exons
awk -v genes="$genes" 'BEGIN{FS=OFS="\t"; n=split(genes, a, " "); for(i=1; i<=n; i++){if(index(a[i], "_introns")!=0){introns[a[i]]=a[i]}else{exon=a[i]"_exons"; exons[exon]=exon}}}{if(introns[$4]==null && exons[$4]==null) next; if(index($4, "introns")!=0){split($4, b, "_"); gene=b[1]; split($3, b, "intron"); name=gene"_intron_"b[3]} else if(index($4, "exons")!=0){split($4, b, "_"); gene=b[1]; split($3, b, "_"gene"_"); name=gene"_exon_"b[2]} print $8, $9, $10, "-", name}' $masterListFile > bedfile.bed.temp


# merge consecutive intervals that are less than 30bp apart
awk 'BEGIN{FS=OFS="\t"; dist=30}{if($1!=chr){if(chr!=null){print chr, start, stop, "-", name;} chr=$1; start=$2; stop=$3; name=$5; next} if($2-stop<dist && $2-stop>=0){ stop=$3;} else if(start-$3<dist && start-$3>=0){start=$2} else{print chr, start, stop, "-", name; chr=$1; start=$2; stop=$3; name=$5;} }END{print chr, start, stop, "-", name;}' bedfile.bed.temp > bedfile.bed


rm bedfile.bed.temp


# add copy number and fingerprinting probes
awk 'BEGIN{FS=OFS="\t"}{if($4=="FingerprintSNPs" || index($4, "_CopyNumber")!=0) print $8, $9, $10, "-", $3 }' $masterListFile >> bedfile.bed


#
