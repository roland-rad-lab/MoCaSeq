#!/bin/bash

##########################################################################################
##
## CheckReferenceFiles.sh
##
## Checks if all mandatory files from the config are found.
##
##########################################################################################

FileList=$@ # takes all variables from the input list

# set binary flag
allfound=yes

for file in $FileList; do
  if [ ! -f $file ]; then
  # check if file is found
  echo "File not found: $file"
  allfound=no
  elif [ ! -s $file ]; then
  # check if file is not empty
  echo "File is empty: $file"
  allfound=no
  fi
done


if [ $allfound == "no" ]; then
  echo "Some mandatory files have issues! Please check the consistency of your genome reference folder or and check for errors during the genome building process"
  exit 1
fi
