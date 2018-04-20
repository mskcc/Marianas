#! /bin/bash

# make title file from a Manifest file.

# $1 - Manifest file
# $2 - Optional pooled normal: FFPEPOOLEDNORMAL or FROZENPOOLEDNORMAL. Used only if the manifest does not have a pooled normal

titleFile=${1/Manifest/title}

printf "Barcode\tPool\tSample_ID\tCollab_ID\tPatient_ID\tClass\tSample_type\tInput_ng\tLibrary_yield\tPool_input\tBait_version\tGender\tPatientName\tMAccession\tExtracted_DNA_Yield\n" > $titleFile

awk 'BEGIN{FS=OFS="\t"}{if($1=="CMO_SAMPLE_ID" || $1=="CMO-SAMPLE-ID") next; sex=$11; if(sex!="Male" && sex!="Female") {sex="-"} print $12, $17, $1, $3, $2, $6, $8, $14, $15, $16, $19, sex, "-", "-", "-"}' $1 >> $titleFile

# add a pooled normal if one does not exist
pooledNormalDir=/ifs/work/bergerm1/Innovation/StandardNormals

count=`cut -f 6 $titleFile | grep -c -e "NormalPool"`

case $count in
"1")
	# we are good, nothing to be done
	;;
"0")
	# add pooled normal
	if [ "$2" == "FFPEPOOLEDNORMAL" ] || [ "$2" == "FROZENPOOLEDNORMAL" ]
	then
		# get the pool value to insert in the pooled normal line
		pool=`tail -n +2 $titleFile | head -1 | cut -f 2`
		# insert the pooled normal line with proper pool value replacement
		tail -n +2 $pooledNormalDir/$2/title_file.txt | awk -v p=$pool 'BEGIN{FS=OFS="\t"}{$2=p; print}' >> $titleFile
	else
		echo -e "pooled normal not supported: $2"
		exit
	fi
	;;
*)
	# throw error
	echo -e "More than 1 pooled normals!"
	exit
	;;
esac


# remove any space characters from the title file
sed -i 's/ //g' $titleFile

# change "NormalPool" to "PoolNormal", only accept the first pooled normal
awk 'BEGIN{FS=OFS="\t"; gotPooledNormal=0;}{if($6=="NormalPool"){if(gotPooledNormal==1){next;} else{$6 = "PoolNormal"; gotPooledNormal=1; print $0;}} else {print $0}}' $titleFile > t
mv t $titleFile


# these checks must be performed at the very end
# check for duplicate sample names
s=`cut -f 3 $titleFile | sort | uniq -d | wc -l`
if [ "$s" != "0" ]
then
	# throw error
	echo -e "Duplicate sample names!!"
	exit
fi

# check for duplicate barcodes
s=`cut -f 1 $titleFile | sort | uniq -d | wc -l`
if [ "$s" != "0" ]
then
	# throw error
	echo -e "Duplicate barcodes!!"
	exit
fi
