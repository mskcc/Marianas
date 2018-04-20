#!/bin/bash

# this script is called from read-file.sh and this action will be perfomred for each line in the file

# $1 - line from the file

#echo "$1"

# split a line on tab
#IFS=$'\t' read -ra array <<< "$1"
#for i in "${array[@]}"
#do
#	# process "$i"
#done


sites=`awk 'BEGIN{FS=OFS="\t"}{if($4<500) next; max=1; count=0; for(i=5; i<=8; i++){if($i>$max) max=i} for(i=5; i<=8; i++){if(i==max) continue; f=$i/$4; if($i>=2 && f >= 0.001 && f < 0.25) count++} if(count>0) print}' "$1" | wc -l`

echo -e "$1\t$sites"













