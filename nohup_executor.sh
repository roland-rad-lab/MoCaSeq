# please make sure nextflow and charliecloud are loaded (see LRZ.md in MoCaSeq repo)
nextflow run -r human-pipeline-nextflow-2 roland-rad-lab/MoCaSeq -profile charliecloud -work-dir /gpfs/scratch/pn29ya/${USER}/${USER}/test/work -entry MAP \
	--input https://raw.githubusercontent.com/roland-rad-lab/MoCaSeq/human-pipeline-nextflow-2/input/remap/human_remap_10p_test.tsv \
	--output_base /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/test/results \
	--genome_build.human GRCh38.p12 \
	--custom_config_version mocaseq-lrz \
	--custom_config_base https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/nextflow-configs \
	--debug true
