#!/bin/bash

# read the contents of the given file line by line and print them

# $1 - the file to be read

while IFS='' read  -r line || [[ -n "$line" ]]
do
	~/software/bin/do-per-line.sh "$line"
done < "$1"


