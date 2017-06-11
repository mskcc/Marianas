#! /bin/bash

# $1 - file with probe pool names, one per line

masterListFile=~/resources/IDT-probes/IDT_CaptureProbe_MASTERLIST_v3.txt

awk -v poolsFile="$1" 'BEGIN{FS=OFS="\t"; while(getline < poolsFile){pools[$1]=$1}}{if(pools[$4]==null) next;  print $8, $9, $10, "-", $3}' $masterListFile > bedfile.bed

# adjust the 4 FP intervals where the snp is actually outside of the interval
awk 'BEGIN{FS=OFS="\t"}{if($5=="FP_rs10102929"){$2=18464920} else if($5=="FP_rs1365740"){$2=60277540} else if($5=="FP_rs10899035"){$2=74401755} else if($5=="FP_rs120434"){$3=110697315} print}' bedfile.bed > bedfile.bed.temp

mv bedfile.bed.temp bedfile.bed


#rm bedfile.bed.temp



#
