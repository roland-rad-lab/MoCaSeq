#!/bin/bash

##########################################################################################
##
## CheckReferenceFiles.sh
##
## Checks if all mandatory files from the config are found.
##
##########################################################################################

FileList=$@ # takes all variables from the input list 

allfound=yes
for file in $FileList; do
  if [ ! -f $file ]; then
  echo "File not found: $file"
  allfound=no
  fi
done

if [ $allfound == "no" ]; then
  echo "Some mandatory files are missing! Please check the consistency of your genome reference folder or and check for errors during the genome building process"
  exit 1
fi
