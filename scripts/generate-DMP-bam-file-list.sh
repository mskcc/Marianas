#! /bin/bash

# A file called PatientID-DMPID.txt must be present in the working directory. It's a tab delimited file where coulmn 1 is patient ID and column 2 is DMP ID

keyFile=/ifs/dmprequest/12-245/key.txt
bamPath=/ifs/dmpshare/share/irb12_245


# Join the PatientID-DMPID and DMPID-BamStem files to get all the necessary bam file stems and corresponding pateint ids that is the input for bam -> standard bam script
awk 'BEGIN{FS=OFS="\t"; while(getline < "PatientID-DMPID.txt"){if(NR==1) continue; split($2, a, ", "); for(i=1; i<=length(a); i++){split(a[i], b, "-"); DMPNumber=b[2]; pid[DMPNumber]=$1;}} FS=","; OFS="\t"}{ split($1, b, "-"); if(b[4]=="IM5"){panel="IMPACT410"} else {panel="IMPACT341"}  if(pid[b[2]]!=null){print $2, pid[b[2]], panel, $1}}'  $keyFile | sort -k 2,2 > t

# for DMP normals, prefer IM5 bams
awk 'BEGIN{FS=OFS="\t"; while(getline < "t"){split($4, a, "-"); if(a[4]=="IM5"){IM5[$2]=$2}}}{ split($4, a, "-"); if(index($1, "-N") != 0 && a[4]=="IM3" && IM5[$2] != null ) next; print }' t > aggregated-info.txt

# create a tab delimited file of full paths to both tumor and normal bams and pateint ids from an aggregated info file as generated above
awk -v path=$bamPath 'BEGIN{FS=OFS="\t";}{ if($1=="#N/A" || seen[$1] != null) next; print path"/"$1".bam", $2, $3; seen[$1] = $1; }' aggregated-info.txt > file-list.txt









#
