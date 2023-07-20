# MoCaSeq launch script for sample group batch02-PCSI_0357
 # path to DSS project folder
 projectDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS
 # path to MoCaSeq git repo
 repoDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq
 # path to working dir
 workDir=/gpfs/scratch/pn29ya/$USER/mocaseq-nextflow/mocaseq
 mkdir -p $workDir
 # nextflow needs this variable to get SLURM job status
 export SLURM\_CLUSTERS="serial"

 nextflow run ${repoDir}/main.nf -profile charliecloud,slurm -work-dir ${workDir} --output_base ${projectDir}/output --genome_build.human GRCh38.p12 --custom_config_version serial-std --custom_config_base ${repoDir}/conf --input batch02-PCSI_0357_Pa_P_5262.tsv
