## Running MoCaSeq nextflow on the LRZ

This aims to get you up and running on the LRZ HPC from scratch. First lets start with some general things that are helpful.

#### Getting started.
There is in fact a lot of software available on the LRZ HPC system.
```bash
# Show available softwares
module avail
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
export PATH="${PATH}:${HOME}/software/charliecloud-0.26/bin"
which ch-run
# /dss/dsshome1/lxc0C/ge26baf2/software/charliecloud-0.26/bin/ch-run (your path will be in your home)


# Nextflow

cd $HOME/software/bin
curl -fsSL https://get.nextflow.io | bash
chmod u+x $HOME/software/bin/nextflow
# $HOME/software/bin should also be on our path
export PATH="${PATH}:${HOME}/software/bin"
which nextflow
# /dss/dsshome1/lxc0C/ge26baf2/software/bin/nextflow (your path will be in your home)
```

#### Where to put things
We have access to various storage systems on the LRZ:
 - Project dir (/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000)
   - 10 Tb
   - This is shared by different users within the same project
   - You are limited in the number of files you can create here so its best for archiving
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

# directory for TCGA project
/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0006/

###### user home dir ######################################
# Extracted container images (You have your own copy so you can debug and hack it without affecting others)
# for my co-workers everyone can use /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images-live instead
$HOME/images-live
# Default location for config files for our pipelines (e.g. MoCaSeq should read $HOME/nextflow-configs/mocaseq/pipeline/mocaseq.config by default)
$HOME/nextflow-configs
# If you wrote a shell script to invoke the pipeline then somewhere in your home dir is a good place for it (e.g. $HOME/pipelines/compass/bin/run.sh)

###### scratch dir #######################################
# /tmp on the LRZ nodes is tiny so instead we use this folder on the scratch space
/gpfs/scratch/pn29ya/${USER}/${USER}/tmp
# folder for nextflow work dir (e.g. compass samples)
/gpfs/scratch/pn29ya/${USER}/${USER}/compass/work
# Technically you could set the results dir to be in the project dir (or copy the results there once the pipeline is complete)

```

### Deploy the pipeline
Nextflow will download and cache the pipeline code directly from the gihub repo, therefore we only need to supply the configuration (tailored for our compute environment) and the software used by the pipeline (we will use container images). We store our [images](https://github.com/roland-rad-lab/Cluster/blob/main/Images.md) on our [LRZ GitLab](https://gitlab.lrz.de/roland-rad-lab/images-public/container_registry), if you want to know more, here are the [docs](https://docs.gitlab.com/ee/user/packages/container_registry/)). Outside of the LRZ we can pull the images we need and create tarballs to copy to LRZ.
```bash

```
Now we can login to LRZ, extract our images and try to configure and run the pipeline.

```bash
# Extracting our images (It is not possible to put them in the project dir as you use up all of your inode quota)
# df -hi /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images-live
# Filesystem     Inodes IUsed IFree IUse% Mounted on
# dssfs02          215K  215K   136  100% /dss/dssfs02

# for me it anyway makes more sense to have your own copy in your home folder so you can experiment without affecting other users
# The tarballs are in /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images

mkdir ${HOME}/images-live
cd ${HOME}/images-live
mkdir mocaseq2 structural-variation-jabba cnv-kit-0.9.9

tar -xzf /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images/mocaseq2.tar.gz -C mocaseq2
tar -xzf /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images/structural-variation-jabba.tar.gz -C structural-variation-jabba
tar -xzf /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images/quay.io%biocontainers%cnvkit:0.9.9--pyhdfd78af_0.tar.gz -C cnv-kit-0.9.9

# The final thing we need to do is to ensure that we have mount points for important resources
# (At some point this should be incorporated into the container images)
mkdir -p mocaseq2/var/pipeline
mkdir -p structural-variation-jabba/var/pipeline
mkdir -p cnv-kit-0.9.9/var/pipeline

# Now we need to pretend there is a HOME for cnvlib
echo -e "HOME=/home/fake\n" >> $HOME/images-live/cnv-kit-0.9.9/ch/environment

# There are example configuration files for LRZ in the test-datasets repo
# you can save them locally or try and use them directly

# On the LRZ the SLURM scheduling system has separate areas called clusters
# for SLURM commands to work we need to specify the area we want to address
# in this case the 'serial' cluster
# more cluster info at https://doku.lrz.de/display/PUBLIC/Available+SLURM+clusters+and+features
# we also set TMPDIR because /tmp on the nodes is too small, even to run
# our container system
export SLURM\_CLUSTERS="serial"
export TMPDIR="/gpfs/scratch/pn29ya/${USER}/${USER}/tmp"
mkdir -p ${TMPDIR}

# If you have not added $HOME/software/bin and $HOME/software/charliecloud-0.26/bin to your path then the following line does that
export PATH="${HOME}/software/bin:${HOME}/software/charliecloud-0.26/bin:${PATH}"


# Here we download a tiny test genome and annotation (tiny.human)
# We also add the --tiny flag to skip steps that break with too little data
# Note that here we only use the charliecloud profile, so nothing will be submitted to slurm (see the next example)
# Please remember to add the slurm profile for real data so you don't start trying to run mutect on the head node

mkdir -p /gpfs/scratch/pn29ya/${USER}/${USER}/test

# this execution of the pipleine does not work
# regarding config files: the pipeline call automatically completes the dir "mocaseq-lrz/pipeline"
# please make sure to have the branch name "mocaseq-nextflow" in the url
# complete file will be like https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/nextflow-configs/mocaseq-lrz/pipeline/mocaseq.config
nextflow run \
	roland-rad-lab/MoCaSeq \
	-r human-pipeline-nextflow-2 \
	-profile charliecloud \
	-work-dir /gpfs/scratch/pn29ya/${USER}/${USER}/test/work \
	--output_base /gpfs/scratch/pn29ya/${USER}/${USER}/test/results \
	--custom_config_version mocaseq-lrz \
	--custom_config_base https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/nextflow-configs \
	--test_config_genome_base https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/nextflow-configs \
	--test_config_genome_version mocaseq \
	--genome_build.human tiny.human \
	--tiny \
	--input https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/testdata/bam/human_design.tsv
	
# after altering the pipeline, please make sure to make a git pull from the lrz where the git repo resides.
# the location can be determined with nextflow info roland-rad-lab/MoCaSeq
# should be something like ${HOME}/.nextflow/assets/roland-rad-lab/MoCaSeq 
```

### Testing with real samples
We can download the Texas Open Cancer Genomes samples from the EGA archive. This will also demonstrate running the [pyega3 client](https://github.com/EGA-archive/ega-download-client) for downloading files from the archive.
```bash
# Although a pre-built image exists, they have not updated it in a long time, so it has (https://github.com/EGA-archive/ega-download-client/issues/139)
# On workstation-05 we can build our own image from git

cd ~/Documents/workspace-docker/charlie-cloud/
./build-image-generic.sh pyega3-fixed Dockerfile.pyega3
scp target/pyega3-fixed.tar.gz ge26baf2@lxlogin2.lrz.de:/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images/

```
Next we connect to LRZ and use the pyega3 client to download some test data.


```bash
# ssh ge26baf2@lxlogin2.lrz.de

mkdir ~/images-live
cd ~/images-live
mkdir pyega3-fixed
tar -xzf /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/images/pyega3-fixed.tar.gz -C pyega3-fixed
mkdir ~/test_open_genomes
cd ~/test_open_genomes
mkdir bin
# We need somewhere to save the data for this project
mkdir -p /dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/TestGenomes/ega


cat > bin/cc.sh <<EOF
#!/usr/bin/env bash
image_name="pyega3-fixed"
data_dir="/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/TestGenomes"

ch-run --no-home -w ${HOME}/images-live/\${image_name} -- bash -c "mkdir -p "\${data_dir}" "$PWD"";ch-run --no-home --unset-env="*" -w --set-env=${HOME}/images-live/\${image_name}/ch/environment --no-passwd --bind "\${data_dir}":"\${data_dir}" --bind "${PWD}":"${PWD}" -c "$PWD" ${HOME}/images-live/\${image_name} -- /bin/bash -c "./bin/ega_download.sh"

EOF

cat > bin/ega_download.sh <<EOF
#!/usr/bin/env bash

output_base="/dss/dssfs02/lwp-dss-0001/pn29ya/pn29ya-dss-0000/projects/TestGenomes/ega"

while IFS= read -r line;
do
	pyega3 -c 30 -cf RadLab_credentials_ega.json fetch \${line} --output-dir \${output_base}
done < egaf_ids.txt

EOF

cat > egaf_ids.txt <<EOF
EGAF00002233848
EGAF00002239689
EOF

cat > RadLab_credentials_ega.json <<EOF
{
    "username": "roland.rad@tum.de",
    "password": "SUPER_SECRET_PASSWORD_HERE"
}
EOF

chmod u+x bin/cc.sh
chmod u+x bin/ega_download.sh

# Now download these files
# this will take some time so run in a screen session
# unfortunately the download service is not so great so often many retries will be needed.
./bin/cc.sh

```

### Remapping to GRCh38.p12
Sometimes the data downloaded are in an older genome build, so the first task is to run the MAP workflow, to remap the BAM files to the current reference.

```bash

# Real samples, this will take a long time so run in a screen session

```
