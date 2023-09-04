# MoCaSeq launch script for sample groupbatch01-RAMP_0007
# new DSS path for project specific long ther storage (save results here)
export DSS2=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS
# git repository path
export REPO=${DSS2}/software/MoCaSeq
workDir=/gpfs/scratch/pn29ya/ga89tog2/mocaseq-nextflow/mocaseq
nextflow run ${REPO}/main.nf -profile charliecloud,slurm -work-dir ${workDir}/mocaseq/work --output_base ${DSS2}/output --genome_build.human GRCh38.p12 --custom_config_version serial-std --custom_config_base ${REPO}/conf --input ${REPO}/input/mocaseq/batch01-RAMP_0007.tsv -resume
