# MoCaSeq launch script for sample groupbatch01-RAMP_0008
 projectDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS
 workDir=/gpfs/scratch/pn29ya/ga89tog2/mocaseq-nextflow
 nextflow run ${projectDir}/MoCaSeq/main.nf -profile charliecloud,slurm -work-dir ${workDir}/mocaseq/work --output_base ${projectDir}/output --genome_build.human GRCh38.p12 --custom_config_version mocaseq-lrz --custom_config_base ${projectDir}/MoCaSeq/conf --input ${projectDir}/MoCaSeq/input/mocaseq/batch01-RAMP_0008/.tsv