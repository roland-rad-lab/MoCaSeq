#!/bin/bash

##########################################################################################
##
## Preparation_GetExemplaryData.sh
##
## Downloads WES and WGS data for animal S821 from https://doi.org/10.1038/nature25459
##
##########################################################################################

download=$1

#md5 WES
# S821-WES.Normal.R1.fastq.gz	2e5fdfeed3054868667af81456abab3c
# S821-WES.Normal.R2.fastq.gz	ec6b00b9b44dd55ec7047f8fedbc498b
# S821-WES.Tumor.R1.fastq.gz	df62e7af80478b65ceb3c3f771fa8eb3
# S821-WES.Tumor.R2.fastq.gz	ddd90a70fa84df84fa82a7f289ca2a89

# md5 WGS
# S821-WGS.Normal.R1.fastq.gz	c1dbb3bb4d2c0e9d1b3c651e9315c1e0
# S821-WGS.Normal.R2.fastq.gz	8fee1ba17688f9ace66978fd1e19bc39
# S821-WGS.Tumor.R1.fastq.gz	1d415615ad844ccc7af32da1d9f81e69
# S821-WGS.Tumor.R2.fastq.gz	0b802662bc092ee48b4bffb70341b745

if [ $download = 'all' ] || [ $download = 'WES' ]; then

	# https://www.ebi.ac.uk/ena/data/view/ERS2066588
	wget -O S821-WES.Normal.R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/006/ERR2230866/ERR2230866_1.fastq.gz
	wget -O S821-WES.Normal.R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/006/ERR2230866/ERR2230866_2.fastq.gz
	# https://www.ebi.ac.uk/ena/data/view/ERS2066550
	wget -O S821-WES.Tumor.R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/008/ERR2230828/ERR2230828_1.fastq.gz
	wget -O S821-WES.Tumor.R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR223/008/ERR2230828/ERR2230828_2.fastq.gz

elif [ $download = 'all' ] || [ $download = 'WGS' ]; then
	# https://www.ebi.ac.uk/ena/data/view/ERS1980539
	wget -O S821-WGS.Normal.R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR221/009/ERR2210079/ERR2210079_1.fastq.gz
	wget -O S821-WGS.Normal.R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR221/009/ERR2210079/ERR2210079_2.fastq.gz
	# https://www.ebi.ac.uk/ena/data/view/ERS1980538
	wget -O S821-WGS.Tumor.R1.fastq.gz  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR221/008/ERR2210078/ERR2210078_1.fastq.gz
	wget -O S821-WGS.Tumor.R2.fastq.gz  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR221/008/ERR2210078/ERR2210078_2.fastq.gz
fi