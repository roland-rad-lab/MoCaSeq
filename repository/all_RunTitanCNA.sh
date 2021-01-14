#!/bin/bash

##########################################################################################
##
## all_RunTitanCNA.sh
##
## Loops over TitanCNA and TitanCNASolution.
##
##########################################################################################

name=$1
repository_dir=$2
threads=$3
sequencing_type=$4
species=$5

numCores=$threads
numClusters=5

cd $name/results/Titan/

if [ $sequencing_type = 'WES' ]; then
	alphaKHigh=2500
	alphaK=2500
elif [ $sequencing_type = 'WGS' ]; then
	alphaKHigh=10000
	alphaK=10000
fi

if [ $species = 'Human' ]; then
	titanchroms="c(1:22)"
elif [ $species = 'Mouse' ]; then
	titanchroms="c(1:19)"
fi

## run TITAN for each ploidy (2,3,4) and clusters (1 to numClusters)
echo "Maximum number of clusters: $numClusters";
for ploidy in $(seq 2 4)
	do
	outDir=run_ploidy$ploidy
	mkdir $outDir
	for numClust in $(seq 1 $numClusters)
	do
	echo "Running TITAN for $numClust clusters.";
	echo "Running for ploidy=$ploidy";
	Rscript $repository_dir/all_TitanCNA.R --id $name --hetFile $name.hetFile.txt --cnFile $name.cnFile.txt \
	--numClusters $numClust --numCores $numCores --normal_0 0.5 --ploidy_0 $ploidy \
	--chrs $titanchroms --estimatePloidy TRUE --outDir $outDir --alphaKHigh=$alphaKHigh --alphaK=$alphaK --species $species
	done
	echo "Completed job for $numClust clusters."
done

Rscript $repository_dir/all_TitanCNASolution.R --ploidyRun2=run_ploidy2 \
--ploidyRun3=run_ploidy3 --ploidyRun4=run_ploidy4 --threshold=0.05 \
--outFile $name.optimalClusters.txt
