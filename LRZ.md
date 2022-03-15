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

```



