# LRZ account prep guide
This document guides you through the steps to prepare an LRZ account for using MoCaSeq nextflow pipeline.
Important paths are
```
# new DSS/project path
/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS
```


## Software installation
Required software are java, nextflow, charliecloud and MoCaSeq human-pipeline-nextflow-2 clone.
You can use a shared version of this software at
```bash
/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS/software
```

To install the software to your $HOME please run the following code
```bash
# create software directory
mkdir -p $HOME/software/bin
cd $HOME/software

# install java (alternatively use more recent java version, but check nextflow compatibilty)
wget https://download.oracle.com/java/19/archive/jdk-19.0.2_linux-x64_bin.tar.gz
tar -xf jdk-19.0.2_linux-x64_bin.tar.gz

# install charlie cloud
wget https://github.com/hpc/charliecloud/releases/download/v0.28/charliecloud-0.28.tar.gz
tar -xzf charliecloud-0.28.tar.gz
cd charliecloud-0.28

# deactivate intel modules for compiling charliecoud, otherwise intel icc will be used and result in errors
module list
# your modules may be differn, please unload all intel modules
module unload intel-mpi/2019-intel intel-mkl/2019 intel/19.0.5

# configure charliecloud compilation
./configure

# compile
make


# install nextflow
cd $HOME/software/bin
curl -fsSL https://get.nextflow.io | bash
chmod u+x $HOME/software/bin/nextflow
```

## Charliecloud containers
To run MoCaSeq we need the exact software for reproducibility.
Therefore we use charliecould containers on LRZ (or Docker on other Platforms)
Run following commands to create the charliecloud containers for MoCaSeq in your $HOME under images-live.
```bash
# make foler for containers
mkdir ${HOME}/images-live
cd ${HOME}/images-live
mkdir mocaseq2 structural-variation-jabba cnv-kit-0.9.9

# extract containers
tar -xzf /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/images/mocaseq2.tar.gz -C mocaseq2
tar -xzf /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/images/structural-variation-jabba.tar.gz -C structural-variation-jabba
tar -xzf /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/images/quay.io%biocontainers%cnvkit:0.9.9--pyhdfd78af_0.tar.gz -C cnv-kit-0.9.9

# create mount points for pipeline
mkdir -p mocaseq2/var/pipeline
mkdir -p structural-variation-jabba/var/pipeline
mkdir -p cnv-kit-0.9.9/var/pipeline

# Now we need to pretend there is a HOME for cnvlib
echo -e "HOME=/home/fake\n" >> $HOME/images-live/cnv-kit-0.9.9/ch/environment
```

## .bashrc
For using the software, having names for important paths and cluster specification, add the following lines to your `.bashrc` file in $HOME.

```bash
# new DSS path for project specific long ther storage (save results here)
export DSS2=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/projects/hPDAC/ICGC_PACA_CA_WGS
# old DSS path from previous MoCaSeq work with results from batch01
export DSS=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS
# scratch path
export SRTCH=/gpfs/scratch/pn29ya/$USER/
# git repository path
export REPO=${DSS2}/software/MoCaSeq/
# directory for temporary files (/gpfs/scratch/ might be cleaned unannounced, don't save results here)
export TMPDIR="/gpfs/scratch/pn29ya/${USER}/tmp"
# path to small debug bam files
export DEBUG_BAM=/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/hPDAC/ICGC_PACA_CA_WGS/ega_download/debug_bams/

# add java, nextflow and charlycloud to PATH
export PATH="${DSS2}/software/jdk-19.0.2/bin/:${PATH}:${DSS2}/software/charliecloud-0.28/bin:${DSS2}/software/bin"
# java specific changes: CLASSPATH and LD_LIBRARY_PATH for java
# export CLASSPATH=${CLASSPATH}:${DSS2}/software/jdk-19.0.2/lib
# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${DSS2}/software/jdk-19.0.2/lib
export JAVA_HOME=${DSS2}/software/jdk-19.0.2
export JAVA_CMD=${DSS2}/software/jdk-19.0.2/bin/java

# the following things are just for conveniance
# change bash apprearance
export PROMPT_DIRTRIM=1

# aliases 
alias ll='ls -la'

# SLURM
# default SLURM cluster one of ("cm2", "cm2_tiny", "serial", "mpp3")
# it is important to set SLURM_CLUSTERS, so nextflow can get slurm status
export SLURM\_CLUSTERS="serial"
alias fairshare='sshare --clusters=serial --user=$USER'
# this command gives you details on a certain job, please provide a job number like: jobdetails xxxxxxx
alias jobdetails='sacct --format=JobID%8,JobName%20,user%8,Partition,NNodes%6,AllocCPUS%9,Time,Start,End,State%11,Elapsed,ExitCode,Reason,MaxRSS,Nodelist%18 --clusters=serial,cm2,cm2_tiny,mpp3 -u ${USER} -j $1'
alias sbash="salloc -pcm2_inter --time=01:00:00 -n 1 srun --pty bash -i"
```

























