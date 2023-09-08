# script to save and encrypt finished mocaseq runs
donor=$1
batch=$2

resultsDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS
mocaseqDir=$resultsDir/output/GRCh38.p12
remapDir=$resultsDir/input/GRCh38.p12

repoDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq

# encryption details
skFile=~/c4ghKeys/COMPASS_c4gh.sec
pkFile=~/c4ghKeys/COMPASS_c4gh.pub
export C4GH_PASSPHRASE='Rad1-COMPASS!' # define the password as envir variable

# determine sample names from mocaseq results dir
samples=$(ls $mocaseqDir | grep $donor)

for samp in $samples
do
	echo $samp

	# 1. move bam and FastQC folders from remap
	mv -v $(find $remapDir -name $samp -type d)/results/* $(find $mocaseqDir -name $samp -type d)/results/

	# 2. move to batch folder
	mv -v $(find $mocaseqDir -name $samp -type d) $mocaseqDir/batch0$batch
	cd $mocaseqDir/batch0$batch

	# 4. encrypt sensitive files
	if [ $(echo $samp | grep '_R') == $samp ]
	then
	echo "Encrypting Normal sample: $samp"
	$repoDir/encrypt_LRZ_results.sh $samp 'Rad1-COMPASS!' "Normal"
	elif [ $(echo $samp | grep -E '_(P|M|X)') == $samp ]
	then
	
	# 3. collect LOH, Mutect and HMMCopy results for purity analysis
	# TODO: how to handle bash Postprocessing?
	
	tar czf ${samp}_PURITYDATA.tar.gz \
	     ${samp}/results/LOH/${samp}.VariantsForLOH.txt \
	     ${samp}/results/Mutect2/${samp}.matched.Mutect2.NoCommonSNPs.txt \
	     ${samp}/results/Mutect2/${samp}.matched.Mutect2.txt \
	     ${samp}/results/Mutect2/${samp}.matched.Mutect2.NoCommonSNPs.OnlyImpact.txt \
	     ${samp}/results/Mutect2/${samp}.matched.Mutect2.NoCommonSNPs.OnlyImpact.CGC.txt \
	     ${samp}/results/HMMCopy/${samp}.HMMCopy.20000.log2RR.txt \
	     ${samp}/results/HMMCopy/${samp}.HMMCopy.20000.segments.txt

	~/.local/bin/crypt4gh encrypt --sk $skFile --recipient_pk $pkFile < ${samp}_PURITYDATA.tar.gz > ${samp}_PURITYDATA.tar.gz.c4gh
	
	# check encryption
	file_to_check=${samp}_PURITYDATA.tar.gz
	if [ -n "$(find "$file_to_check.c4gh" -prune -size +1c)" ]; then
        	echo "File successfully encrypted: ${file_to_check}"
        	rm ${file_to_check}
    	else
        	echo "File was NOT successfully encrypted: ${file_to_check}"
        	echo "This file was NOT encrypted! (maybe the password in argument 2 is wrong?)"
        	rm ${file_to_check}.c4gh # remove empty file 
    	fi
	
	echo "Encrypting Tumor sample: $samp"
	$repoDir/encrypt_LRZ_results.sh $samp 'Rad1-COMPASS!' "Tumor"
	else
	echo "Error: cannot detect sample type in $samp. Will not encrypt."
	fi
done
