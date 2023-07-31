# specify output/project dir
projectDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS

# define MoCaSeq version
mocaseqDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq

# specify working dir
workingDir=/gpfs/scratch/pn29ya/$USER/mocaseq-nextflow/remap


# nextflow call with charliecloud profile 
nextflow run ${mocaseqDir}/main.nf -profile charliecloud -entry MAP -work-dir ${workingDir} --output_base ${projectDir}/input --genome_build.human GRCh38.p12 --custom_config_base ${mocaseqDir}/conf --custom_config_version serial-std --input batch02-PCSI_0683.tsv -resume

