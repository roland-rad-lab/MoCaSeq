# MoCaSeq launch script for sample group batch04-RAMP_0005
projectDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS
workDir=/gpfs/scratch/pn29ya/$USER/mocaseq-nextflow/remap
nextflow run /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/main.nf -profile charliecloud,slurm -entry MAP -with-report -with-timeline -work-dir ${workDir}/remap/work --output_base ${projectDir}/input --genome_build.human GRCh38.p12 --custom_config_version serial-std --custom_config_base /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/conf --input batch01-RAMP_0007.tsv
