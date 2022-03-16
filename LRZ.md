## Running MoCaSeq nextflow on the LRZ

This aims to get you up and running on the LRZ HPC from scratch. First lets start with some general things that are helpful.

#### Getting started.
There is in fact a lot of software available on the LRZ HPC system.
```bash
# Show available softwares
moduel avail
# See what versions of specific software (nextflow and charliecloud) are offered
module avail nextflow charliecloud

# To add the sofware to your path you could do
module load charliecloud/0.22

```

To ask for these softwares to be upgraded or for other problems you can sumit a ticket to the LRZ [helpdesk](https://servicedesk.lrz.de). Right now we would like to have at least:
 - charliecloud 0.26 (Because of [#1117](https://github.com/hpc/charliecloud/pull/1117) the tiny TMPDIR on LRZ nodes would cause your jobs to crash)
 - nextflow 21.10.0 (Because we use DSL2 and its fairly new so getting some bug fixes)

At the moment we should install our own version of these two softwares, pending a request to the helpdesk to make some updated version available. This can be done in our home directory, as follows:

```bash

# https://hpc.github.io/charliecloud/install.html
mkdir -p $HOME/software/bin
cd $HOME/software
curl -L https://github.com/hpc/charliecloud/releases/download/v0.26/charliecloud-0.26.tar.gz > charliecloud-0.26.tar.gz
tar -xzf charliecloud-0.26.tar.gz
cd charliecloud-0.26

# Charliecloud will barf if it finds you are trying to use the intel icc as our compiler. On LRZ you are probably using this by default
# We can see what we have loaded using:
module list
# Currently Loaded Modulefiles:
# 1) admin/1.0   2) tempdir/1.0   3) lrz/1.0   4) spack/21.1.1   5) intel/19.0.5   6) intel-mkl/2019   7) intel-mpi/2019-intel  

# So we need to get rid of the intel stuff
module rm intel-mpi/2019-intel intel-mkl/2019 intel/19.0.5 

./configure 
# All we care about here is that we can compile ch-run
# Test suite
# ~~~~~~~~~
# 
#   basic tests, all stages: no
#     test suite enabled ... yes
#     ch-run(1) ... yes

make

# Now we should add $HOME/software/charliecloud-0.26/bin to our path
which ch-run
# /dss/dsshome1/lxc0C/ge26baf2/software/charliecloud-0.26/bin/ch-run


# Nextflow

cd $HOME/software/bin
curl -fsSL https://get.nextflow.io | bash
# $HOME/software/bin should also be on our path
which nextflow
# /dss/dsshome1/lxc0C/ge26baf2/software/bin/nextflow
```

#### Where to put things
We have access to various storage systems on the LRZ:
 - Project dir (/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000)
   - 10 Tb
   - This is shared by different users within the same project
 - Scratch (/gpfs/scratch/pn29ya/${USERNAME})
   - Space on a 2 Pb partition (unlimited)
   - Old files are automatically deleted (so don't put scripts here, and save results you want to keep to the Project dir)
 - Home dir (${HOME})
   - Space on a 0.5 Pb partition (unlimited)

I have the following folders setup to run the MoCaSeq nextflow pipeline:
```bash
####### project dir #########################################
# Container image tarballs
/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images
# MoCaSeq reference folders and other genome references
/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/reference
# Results and input files (safe from deletion)
/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects

###### user home dir ######################################
# Extracted container images (You have your own copy so you can debug and hack it without affecting others)
$HOME/images-live
# Default location for config files for our pipelines (e.g. MoCaSeq should read $HOME/nextflow-configs/mocaseq/pipeline/mocaseq.config by default)
$HOME/nextflow-configs
# If you wrote a shell script to invoke the pipeline then somewhere in your home dir is a good place for it (e.g. $HOME/pipelines/compass/bin/run.sh)

###### scratch dir #######################################
# /tmp on the LRZ nodes is tiny so instead we use this folder on the scratch space
/gpfs/scratch/pn29ya/ge26baf2/tmp
# folder for nextflow work dir (e.g. compass samples)
/gpfs/scratch/pn29ya/ge26baf2/compass/work
# Technically you could set the results dir to be in the project dir (or copy the results there once the pipeline is complete)

```

### Deploy the pipeline
Nextflow will download and cache the pipeline code directly from the gihub repo, therefore we only need to supply the configuration (tailored for our compute environment) and the software used by the pipeline (we will use container images). We store our [images](https://github.com/roland-rad-lab/Cluster/blob/main/Images.md) on our [LRZ GitLab](https://gitlab.lrz.de/roland-rad-lab/images-public/container_registry), if you want to know more, here are the [docs](https://docs.gitlab.com/ee/user/packages/container_registry/)). Outside of the LRZ we can pull the images we need and create tarballs to copy to LRZ.
```bash

```
Now we can login to LRZ, extract our images and try to configure and run the pipeline.

```bash
# There are example configuration files for LRZ in the test-datasets repo
# you can save them locally or try and use them directly

# On the LRZ the SLURM scheduling system has separate areas called clusters
# for SLURM commands to work we need to specify the area we want to address
# in this case the 'serial' cluster
# we also set TMPDIR because /tmp on the nodes is too small, even to run
# our container system
export SLURM_CLUSTERS="serial"
export TMPDIR="/gpfs/scratch/pn29ya/ge26baf2/ge26baf2"

nextflow run \
	roland-rad-lab/MoCaSeq \
	-r human-pipeline-nextflow \
	-resume \
	-dump-channels \
	-profile charliecloud \
	-work-dir /gpfs/scratch/pn29ya/ge26baf2/ge26baf2/test/work \
	--output_base /gpfs/scratch/pn29ya/ge26baf2/ge26baf2/test/results \
	--custom_config_version mocaseq-lrz \
	--custom_config_base https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/nextflow-configs \
	--input https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/testdata/bam/human_design.tsv


```

