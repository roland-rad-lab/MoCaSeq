[![DOI](https://zenodo.org/badge/162331694.svg)](https://zenodo.org/badge/latestdoi/162331694)
[![Docker](https://img.shields.io/badge/Docker-avaible-green.svg)](https://cloud.docker.com/repository/docker/rolandradlab/natureprotocols2018/general)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Computational pipelines for the analysis of mouse cancers

### Contents

* [Overview](#Overview)
    + [Abstract](#abstract)
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
Mouse models of human cancer have proven invaluable in linking genetics to mechanisms and phenotypes. The tremendous opportunities for reverse as well as forward genetics in mice gain further momentum through the sequencing revolution. Sequencing analysis pipelines were however developed for humans, and hence do not account for species-specific differences in genome structures and experimental setups.

Here, we describe standardised computational pipelines tailored specifically for the analysis of mouse genomic data. We present workflows for the detection of all genetic alterations, including single nucleotide variants, indels, copy number variation, loss of heterozygosity and complex rearrangements. All components have been extensively validated and cross-compared using multiple methodologies. The protocol also contributes novel analytical tools, such as pipelines for inference of chromothripsis.

We provide code for all workflows, give step by step guidance on the conduction of individual analysis types and provide advice for data interpretation. The protocol takes 2-7 days, depending on the desired analyses.

## System Requirements
- Using the *bash* version: Linux, we run this pipeline under Ubuntu 18.04 LTS.
- Using the *Docker* version: Platform of your choice.
- Hardware: 
	- Minimum: 8-core processor, 16 GB RAM
	- Optimal (running multiple samples in parallel): 48-core processor, 256 GB RAM, Solid-State Drive
- Disk space: 
	- Reference files: 15 GB
	- Results: 30 GB (WES 100x), 300 GB (WGS)
	- Temporary files during analysis: ~170 GB (WES), ~1000 GB (WGS)

## Installation
1. This pipeline requires a number of tools used for the detection of somatic mutations. Installation procedures are listed in `repository/Preparation_SystemSetup.sh`.

2. Edit `config.sh` to update the relevant paths for each tool, as well as the temporary directory and the number of CPU threads and RAM to be used.

3. Run `sh repository/Preparation_GetReferenceData.sh config.sh` to automatically download all reference files needed to the current directory (inside a newly created folder `Genomes`).

## Usage
### Inputs and their formats
fastq
bam

### Using bash

The complete pipeline is wrapped inside a shell-script. During the analysis, the fastq-files are automatically copied to the target directory.
```
sh DNA/CancerGenomeAnalysis.sh \
-nf '/media/rad/HDD2/examples/S821/S821-WES.Normal.R1.fastq.gz' \
-nr '/media/rad/HDD2/examples/S821/S821-WES.Normal.R2.fastq.gz'  \
-tf '/media/rad/HDD2/examples/S821/S821-WES.Tumor.R1.fastq.gz'  \
-tr '/media/rad/HDD2/examples/S821/S821-WES.Tumor.R2.fastq.gz' \
--name S821-WES --species Mouse --sequencing_type WES --config /media/rad/SSD1/DNA/configadapted.sh \
--artefact none --filtering all --phred phred33 --Delly no --Chromothripsis no --Mutect2 yes --Titan yes 

sh MouseCancerGenomeAnalysis.sh <name> \
<fastq_normal_1> <fastq_normal_2> <fastq_tumor_1> <fastq_tumor_2> \
<sequencing_type> <config_file> [none|CT|GT] [phred33|phred64]"
```

### Using Docker

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