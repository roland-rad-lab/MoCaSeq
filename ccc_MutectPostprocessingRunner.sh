#!/bin/bash
#SBATCH -o ./debug/%x.%j.%N.out
#SBATCH -e ./debug/%x.%j.%N.err
#SBATCH -D ./
#SBATCH -J MutectPostprocessMonitor
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --reservation=gen_seq

module load slurm_setup

#==========================================================================
# MONITORING
# -> create a log file for each node to be monitored
# -> periodically running "ps" on all compute nodes of this job to get CPU load and memory usage
# -> writing "ps" output (snapshot) to log files for each node
#==========================================================================
NODEMONITOR=1             # enable (1) or disable (0) monitoring
MONPATH=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/log  # storage location for log files
DELAY=60                  # intervall (in s) to take snapshots of table of running processes
MAXSNAPSHOTS=1000         # MAXIMUM NUMBER OF SNAPSHOTS TO AVOID CREATING A HUGE LOG FILE!!!
COMMAND=ch-run             # name of command (running on the node) to be checked

if (( NODEMONITOR != 0 )); then
    # prepare monitoring
    if [ -d $MONPATH ]; then rm -rf $MONPATH/*; else mkdir -p $MONPATH; fi
    mpiexec hostname | sort -u > $MONPATH/mpi_hostfile
    # run monitoring script in the background
    mpiexec -f $MONPATH/mpi_hostfile -ppn 1 ./monitor_ps.sh \
                                              $MONPATH $COMMAND $DELAY $MAXSNAPSHOTS &
fi

#==========================================================================
# RUN USER APPLICATION
#==========================================================================
# EDIT HERE: put all your job script stuff here

# run charliecloud container to do Mutect2 post processing
# set container path
cccDir=${HOME}/images-live/mocaseq2

# specify mount paths
workingDir=/gpfs/scratch/pn29ya/$USER/mocaseq-slurm/mocaseq/
mocaseqDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software/MoCaSeq/
outDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/output/GRCh38.p12/
referencesDir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/references/bashMoCaSeq/
# make sure the working dir exists and create a symbolic link to the references named "ref" inside working directory
mkdir -p $workingDir
mkdir -p /gpfs/scratch/pn29ya/$USER/tmp
ln -s $referencesDir ref
mv ref $workingDir

# postprocessing variables
name=PCSI_0410_Ag_M_526

# MoCaSeq call inside charliecloud container
ch-run $cccDir --no-home --set-env=name=${name} -w --no-passwd \
--bind ${workingDir}:/var/pipeline/ \
--bind ${mocaseqDir}:/opt/MoCaSeq/ \
--bind ${outDir} \
--bind ${referencesDir} \
-c $outDir/batch02 \
-- /bin/bash /opt/MoCaSeq/launch/ccc_Mutect_Postprocessing.sh \
	$name $species $config_file $filtering $artefact_type $GATK $type

#==========================================================================
# TERMINATE MONITORING (IF STILL RUNNING)
#==========================================================================
if (( NODEMONITOR != 0 )); then
    mpiexec -f $MONPATH/mpi_hostfile -ppn 1 killall -u $USER -e monitor_ps.sh
fi

