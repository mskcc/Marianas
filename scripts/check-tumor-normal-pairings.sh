#! /bin/bash

# flag the pairs in MergedPool_NormalUsedInMutationCalling.txt that do not have same patient ID's.
# also flag those pairs where the tumor sample is not actually a tumor or normal sample is not actually a normal.

# read each line from the title file, split on tab (bash is pathetic!)
while IFS=$'\t' read -r -a array
do
	tumorSample=${array[0]}
	normalSample=${array[1]}

	tumor=`grep -w $tumorSample title_file.txt | cut -f 6`
	normal=`grep -w $normalSample title_file.txt | cut -f 6`
	tumorPatientID=`grep -w $tumorSample title_file.txt | cut -f 5`
	normalPatientID=`grep -w $normalSample title_file.txt | cut -f 5`

	if [ "$tumor" != "Tumor" ]
	then
		printf "Not a tumor: $tumorSample\n"
	fi

	if [ "$normal" != "Normal" ]
	then
		printf "Not a normal: $normalSample\n"
	fi

	if [ "$tumorPatientID" != "$normalPatientID" ]
	then
		printf "$tumorSample\t$tumorPatientID\t$normalSample\t$normalPatientID\n"
	else
		printf "Patient IDs match\n"
	fi

done < MergedPool_NormalUsedInMutationCalling.txt