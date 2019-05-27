[![DOI](https://zenodo.org/badge/162331694.svg)](https://zenodo.org/badge/latestdoi/162331694)
[![Docker](https://img.shields.io/badge/Docker-avaible-green.svg)](https://cloud.docker.com/repository/docker/rolandradlab/natureprotocols2018/general)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## MoCaSeq: Computational pipelines for the analysis of mouse cancers

__Sebastian Lange<sup>1,2,3</sup>, Thomas Engleitner<sup>1,2,3</sup>, Sebastian Mueller<sup>1,2,3</sup>, Roman Maresch<sup>1,2,3</sup>, Maximilian Zwiebel<sup>1,2,3</sup>, Laura Gonzalez-Silva<sup>4</sup>, Günter Schneider<sup>2</sup>, George S. Vassiliou<sup>5,6,7</sup>, Mathias J. Friedrich<sup>1,2,3</sup>, Dieter Saur<sup>2,3,8</sup>, Ignacio Varela<sup>4</sup>, Roland Rad<sup>1,2,3,8</sup>__
<br>
<sub>
	1 Institute of Molecular Oncology and Functional Genomics, TUM School of Medicine, Technische Universität München, 81675 Munich, Germany<br>
	2 Department of Medicine II, Klinikum rechts der Isar, Technische Universität München, 81675 Munich, Germany<br>
	3 Center for Translational Cancer Research (TranslaTUM), TUM School of Medicine, Technische Universität München, 81675 Munich, Germany<br>
	4 Instituto de Biomedicina y Biotecnología de Cantabria. Universidad de Cantabria-CSIC. 39012 Santander, Spain<br>
	5 The Wellcome Trust Sanger Institute, Genome Campus, Hinxton, CB10 1SA Cambridge, UK<br>
	6 Wellcome Trust-MRC Stem Cell Institute, Cambridge Biomedical Campus, University of Cambridge, CB2 0SL Cambridge, UK<br>
	7 Department of Haematology, Cambridge University Hospitals NHS Trust, CB2 0QQ Cambridge, UK<br>
	8 German Cancer Consortium (DKTK), German Cancer Research Center (DKFZ), 69120 Heidelberg, Germany<br>
</sub>

### Contents

* [Overview](#Overview)
    + [Abstract](#abstract)
    + [TL;DR](#tl;dr)
* [System Requirements](#system-requirements)
* [Installation](#installation)
* [Usage](#usage)
    - [Inputs and their formats](#inputs-and-their-formats)
    - [Using bash](using-bash)
    - [Using Docker](using-docker)
    - [Example](#example)
* [Bug reports](#bug-reports)
* [Citation](#citation)
* [License](./LICENSE)

## Overview
This repository serves as a companion to an analysis workflow protocol for mouse cancer next-generation data.

### Abstract
Mouse models of human cancer have proven invaluable in linking genetics to molecular mechanisms and phenotypes. The tremendous opportunities for reverse as well as forward genetics in mice gain further momentum through the sequencing revolution. Sequencing analysis pipelines were however developed for humans, and hence do not account for species-specific differences in genome structures and experimental setups. 

Here, we describe standardised computational pipelines tailored specifically for the analysis of mouse genomic data. We present workflows for the detection of all genetic alterations, including single nucleotide variants, indels, copy number variation, loss of heterozygosity and complex rearrangements. All components have been extensively validated and cross-compared using multiple methodologies. The protocol also contributes novel analytical tools, such as pipelines for inference of chromothripsis. 

We give step by step guidance on the conduction of individual analysis types and provide advice for data interpretation. The complete code is available [online](https://github.com/roland-rad-lab/MoCaSeq). 

This protocol takes 2-7 days, depending on the desired analyses.  

## TL;DR
1. Download and install [Docker](https://www.docker.com/products/docker-desktop).
2. Create a new folder used for testing. Use a volume with at least 250 GB of free disk space.
```
mkdir MoCaSeq && cd MoCaSeq
```
2. Replace `<your_working_directory>` and `<mocaseq_version>` run 
```
sudo docker run \
-v <your_working_directory>/:/var/pipeline/ \
rolandradlab/mocaseq:<mocaseq_version> \
--test yes
``` 
This automatically downloads the container, downloads all required reference files and tests the pipeline. If succesful, reference files will be located in `ref/` and test results in `MoCaSeq_Test/`. This will take up to 24 hours!

3. Download and unzip the github repository: 
```
wget https://github.com/rolandradlab/MoCaSeq/archive/master.zip \
&& unzip master.zip \
&& rm master.zip \
&& mv MoCaSeq-master MoCaSeq
```
4. Use the provided script to download both tumor and matched normal fastq-files from one pancreatic cancer, which developed in a conditionally-activated Kras<sup>G12D</sup>-model.
```
sh MoCaSeq/repository/Preparation_GetExemplaryData.sh WES
```
`all` will download both WES (100x) and WGS (30x) data, using 100 GB of disk space. Use `WES` or `WGS` to only download the respective files.

The raw data is available from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) using the run accessions ERR2230828 (WES Tumour), ERR2230866 (WES Normal), ERR2210078 (WGS Tumour) and ERR2210079 (WGS Normal).

5. Replace `<your_working_directory>` and `<mocaseq_version>`. This directory should contain `ref/`, which was generated in step 2. Additionally, replace `<threads>` and `<ram>`, then run the pipeline using the downloaded data Depending on the available CPU and RAM, this will take approximately 24 hours.
```
sudo docker run \
-e USERID=`id -u` -e GRPID=`id -g` \
-v <your_working_directory>:/var/pipeline/ \
rolandradlab/mocaseq:<mocaseq_version>  \
-nf '/var/pipeline/S821-WES.Normal.R1.fastq.gz' \
-nr '/var/pipeline/S821-WES.Normal.R2.fastq.gz' \
-tf '/var/pipeline/S821-WES.Tumor.R1.fastq.gz' \
-tr '/var/pipeline/S821-WES.Tumor.R2.fastq.gz' \
--temp_dir /var/pipeline/temp --config /opt/MoCaSeq/config_docker.sh \
--name S821-WES --species Mouse --threads <threads> --ram <ram> --sequencing_type WES \
--artefact GT --filtering hard --Delly no --Mutect2 yes
```
6. Browse the results located in `S821/results/`.

## System Requirements
- Using the *bash* version: Linux, we run this pipeline under Ubuntu 18.04 LTS.
- Using the *Docker* version: Platform of your choice.
- Hardware: 
	- Minimum: 8-core processor, 32 GB RAM
	- Optimal (running multiple samples in parallel): 48-core processor, 256 GB RAM, Solid-State Drive
- Disk space: 
	- Reference files: 15 GB
	- Results: 30 GB (WES 100x), 300 GB (WGS)
	- Temporary files during analysis: ~170 GB (WES), ~1000 GB (WGS)

## Installation
All *bash*, *R* and *python* scripts are directly invokable from `repository`.

1. This pipeline requires a number of tools used for the detection of somatic mutations. Installation procedures are listed in `repository/Preparation_SystemSetup.sh`.

2. Edit `config.sh` to update the relevant paths for each tool, as well as the temporary directory and the number of CPU threads and RAM to be used.

3. Run `sh repository/Preparation_GetReferenceData.sh config.sh` to automatically download all reference files needed to the current directory (inside a newly created folder `Genomes`).

## Usage
### Inputs and their formats
fastq
bam

### Using Docker

### Docker
We provide a docker image containing the complete analysis pipeline in order to simplify deployment and to keep software versions as consistent as possible. The basic commandline to run the dockerized pipeline is as follows:
```
sudo docker run \
-v <your_working_directory>/:/var/pipeline/ \
rolandradlab/mocaseq:<github_release> <options>
```
Options are identical to calling `CancerGenomeAnalysis.sh` directly and are listed [below](#options). When invoked without options, the container will start, display usage information and shut down.

Docker by design runs the container and its contents as user root (UID 1 and GID 1). Persistent directories mounted into the container with the option '-v' therefore are owned by root. If this is undesirable, you can pass the UID and GID of your current user into the container by specifying ``-e USERID=`id -u` -e GRPID=`id -g` `` when calling `docker run`:
```
sudo docker run \
-v <your_working_directory>/:/var/pipeline/ \
-e USERID=`id -u` -e GRPID=`id -g` \
rolandradlab/mocaseq:<github_release> <options>
```
This will run the pipeline with your current UID and GID and set the permissions of the output files accordingly.

### Using bash

The complete pipeline is wrapped inside a shell-script. During the analysis, the fastq-files are automatically copied to the target directory.
```
sh MoCaSeq/CancerGenomeAnalysis.sh \
-nf '/media/rad/HDD2/examples/S821/S821-WES.Normal.R1.fastq.gz' \
-nr '/media/rad/HDD2/examples/S821/S821-WES.Normal.R2.fastq.gz'  \
-tf '/media/rad/HDD2/examples/S821/S821-WES.Tumor.R1.fastq.gz'  \
-tr '/media/rad/HDD2/examples/S821/S821-WES.Tumor.R2.fastq.gz' \
--name S821-WES --species Mouse --sequencing_type WES --config /media/rad/SSD1/MoCaSeq/configadapted.sh \
--artefact none --filtering all --phred phred33 --Delly no --Chromothripsis no --Mutect2 yes --Titan yes 

sh MouseCancerGenomeAnalysis.sh <name> \
<fastq_normal_1> <fastq_normal_2> <fastq_tumor_1> <fastq_tumor_2> \
<sequencing_type> <config_file> [none|CT|GT] [phred33|phred64]"
```

### Options


### Example
1. Download examplary data to the current folder.
Choose from `WES` (download WES example), `WGS` (download WGS example) or `all` (both WES and WGS).
```
sh $repository_dir/Preparation_GetExemplaryData.sh WES
```

2. Run the main workflow using the downloaded files.

```
sh MouseCancerGenomeAnalysis.sh S821-WES \
S821-WES.Normal.R1.fastq.gz S821-WES.Normal.R2.fastq.gz \
S821-WES.Tumor.R1.fastq.gz S821-WES.Tumor.R2.fastq.gz \
WES config.sh none phred33
```

3. The results can be found in `S821-WES/results/`

## Bug reports
Please send comments and bug reports to: sebastian.lange [@] tum.de

## Citation
Please cite our primary paper: \
Mueller, S., Engleitner, T., Maresch, R., Zukowska, M., Lange, S., …, Rad, R. (2018). \
Evolutionary routes and KRAS dosage define pancreatic cancer phenotypes. Nature, 554(7690), 62–68. \
https://doi.org/10.1038/nature25459

## License
The program is distributed under the MIT license.