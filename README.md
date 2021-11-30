[![DOI](https://zenodo.org/badge/162331694.svg)](https://zenodo.org/badge/latestdoi/162331694)
[![Docker](https://img.shields.io/badge/Docker-avaible-green.svg)](https://hub.docker.com/r/rolandradlab/mocaseq)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## MoCaSeq: Analysis pipelines for cancer genome sequencing in mice

__Sebastian Lange<sup>1,2,3</sup>, Thomas Engleitner<sup>1,3</sup>, Sebastian Mueller<sup>1,3</sup>, Roman Maresch<sup>1,3</sup>, Maximilian Zwiebel<sup>1,3</sup>, Laura Gonzalez-Silva<sup>4</sup>, Günter Schneider<sup>2</sup>, Ruby Banerjee<sup>5</sup>, Fengtang Yang<sup>5</sup>, George S. Vassiliou<sup>5,6,7</sup>, Mathias J. Friedrich<sup>1,2,3</sup>, Dieter Saur<sup>2,3,8,9</sup>, Ignacio Varela<sup>4</sup>, Roland Rad<sup>1,2,3,8</sup>__
<br>
<sub>
	1 Institute of Molecular Oncology and Functional Genomics, School of Medicine, Technische Universität München, 81675 Munich, Germany<br>
	2 Department of Medicine II, Klinikum rechts der Isar, Technische Universität München, 81675 Munich, Germany<br>
	3 Center for Translational Cancer Research (TranslaTUM), School of Medicine, Technische Universität München, 81675 Munich, Germany<br>
	4 Instituto de Biomedicina y Biotecnología de Cantabria. Universidad de Cantabria-CSIC. 39012 Santander, Spain<br>
	5 The Wellcome Trust Sanger Institute, Genome Campus, Hinxton, CB10 1SA Cambridge, UK<br>
	6 Wellcome Trust-MRC Stem Cell Institute, Cambridge Biomedical Campus, University of Cambridge, CB2 0SL Cambridge, UK<br>
	7 Department of Haematology, Cambridge University Hospitals NHS Trust, CB2 0QQ Cambridge, UK<br>
	8 German Cancer Consortium (DKTK), German Cancer Research Center (DKFZ), 69120 Heidelberg, Germany<br>
	9 Institute for Translational Cancer Research and Institute for Experimental Cancer Therapy, School of Medicine, Technische Universität München, 81675 Munich, Germany<br>
</sub>

### Contents

* [Overview](#Overview)
    + [Abstract](#abstract)
* [System Requirements](#system-requirements)
* [Input formats](#input-formats)
* [Usage](#usage)
    - [Available Analyses](#available-analyses)
    - [Input TSV](#input-tsv)
    - [Options](#options)
    - [Folder locations](#folder-locations)
    - [Interactive mode](#interactive-mode)
* [TL;DR](#tldr)
* [Bug reports](#bug-reports)
* [Citation](#citation)
* [License](./LICENSE)

## Overview
This repository serves as a companion to an analysis workflow protocol for mouse cancer next-generation data. The manuscript, published in *Nature Protocols*, can be accessed [online](https://t.co/bRmsqg2lgb?amp=1). The pipeline has been extended to also run human cancer data and now uses [nextflow](https://www.nextflow.io/) for ease of use and increased support of high performance computing environments.

### Abstract
Mouse models of human cancer have transformed our ability to link genetics, molecular mechanisms and phenotypes. Both reverse and forward genetics in mice are currently gaining momentum through advances in next generation sequencing. Methodologies to analyse sequencing data were however developed for humans, and hence do not account for species-specific differences in genome structures and experimental setups.

Here, we describe standardised computational pipelines tailored specifically for the analysis of mouse genomic data. We present novel tools and workflows for the detection of different alteration types, including single nucleotide variants, indels, copy number variation, loss of heterozygosity and complex rearrangements, such as chromothripsis.

Workflows have been extensively validated and cross-compared using multiple methodologies. We also give step by step guidance on the execution of individual analysis types and provide advice on data interpretation.

The complete code is available [online](https://github.com/roland-rad-lab/MoCaSeq/tree/human-pipeline-nextflow) and can be run as a nextflow pipeline.

```bash
# input is a required parameter and is a tab separated file providing information about samples and file paths
nextflow run roland-rad-lab/MoCaSeq -r human-pipeline-nextflow --input input.tsv
```

This protocol takes 2-7 days, depending on the desired analyses.  

## System Requirements
We **strongly** recommend to use the container images that accompany this pipeline and encapsulate the required software.

- This pipeline was developed under under Ubuntu 20.04 LTS.
- To use the provided container images we recommend using [Charliecloud](https://hpc.github.io/charliecloud/index.html) a contiainer runtime that does not require administrator privieges.

- Hardware:
	- Minimum: 8-core processor, 32 GB RAM
	- Optimal (running multiple samples in parallel): 48-core processor, 256 GB RAM, Solid-State Drive
- Disk space:
	- Reference files: 15 GB
	- Results: 30 GB (WES 100x), 300 GB (WGS)
	- Temporary files during analysis: ~170 GB (WES), ~1000 GB (WGS)

Most of the code is written as nextflow `modules` but some is imported as libraries from `repository`.
All *bash*, *R* and *python* scripts are directly invokable from `repository`.

## Input formats
The standard input format are FASTQ files produced from modern Illumina sequencers. These should be split into forward and reverse reads for both the tumour and the matched normal sample. BAM files can be processed as well, giving the user the choice of starting after alignment or re-mapping all the raw data (see below).

## Usage
The pipeline requires a tsv file containing the sample information and file paths to the fastq or bam files. We also recommend supplying a custom configuration with details of your computational environment such as genome reference file locations, scheduling system, resource limits etc.

We provide container images containing the complete software used by the analysis pipeline in order to simplify deployment and to keep software versions as consistent as possible. You can find example configuration files in `example` and in our [test data](https://github.com/roland-rad-lab/test-datasets/tree/mocaseq-nextflow). The configuration is separated into three files:
- mocaseq.config: Containers, resources requirements, task scheduling
- genomes.config: Location of reference genomes and associated index files
- genome\_annotations.config: Location of annotations on the reference genome, e.g. genes, common variants 

Once you have setup mocaseq.config for your computing environment you can use [genomes.config](https://github.com/roland-rad-lab/test-datasets/blob/mocaseq-nextflow/nextflow-configs/mocaseq/pipeline/genomes.config) and [genome\_annotations.config](https://github.com/roland-rad-lab/test-datasets/blob/mocaseq-nextflow/nextflow-configs/mocaseq/pipeline/genome_annotations.config) directly to run the test data.

```bash
# test_config_genome_* can be used with custom_config_* to specify a different location for the genome config files
nextflow run \
	roland-rad-lab/MoCaSeq \
	-r human-pipeline-nextflow \
	-entry MAP \
	--test_config_genome_base https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/nextflow-configs \
	--test_config_genome_version mocaseq \
	--genome_build.human tiny.human \
	--input https://raw.githubusercontent.com/roland-rad-lab/test-datasets/mocaseq-nextflow/testdata/bam/human_remap_design.tsv

```

The following example will look for all three config files in ${HOME}/nextflow-configs/human-pipeline-nextflow/pipeline.


```bash
# To run with the charliecloud profile defined in the custom configuration files in the directory ${HOME}/nextflow-configs/human-pipeline-nextflow you can start the pipeline as follows:
nextflow run \
	roland-rad-lab/MoCaSeq \
	-r human-pipeline-nextflow \
	-profile charliecloud \
	--custom_config_version human-pipeline-nextflow \
	--custom_config_base ${HOME}/nextflow-configs \
	--input input.full.tsv
```

### Available Analyses
Although most analysis is specified in the input file (a mix of mouse and human samples is fine) it is also possible to specify an alternative workflow using the nextflow entry option.
- PON
- REMAP

```bash
# To run the REMAP workflow
nextflow run roland-rad-lab/MoCaSeq -r human-pipeline-nextflow -entry REMAP --input input.tsv
```

### Input TSV
| Column header    | Details |
|------------------|---------|
| Sample_Name      | Unique name for each biological isolate (Must be different for a Tumor/Normal pair) |
| Sample_Group     | Unique name for a group of samples that should be analysed together (i.e. shared by a Tumor/Normal pair) |
| Library_ID       | Sequencing library id (ignored for now) |
| Lane             | Sequencing lane (ignored for now) |
| Colour_Chemistry | Property of sequence data (ignored for now)|
| SeqType          | 'wgs' for genome or 'wex' for exome data |
| Organism         | 'human' or 'mouse' |
| Type             | 'Normal' or 'Tumor' |
| R1               | FASTQ file for read 1 (or BAM file path if remapping) |
| R2               | FASTQ file for read 2 if paired end (or repeat BAM file path if remapping paired data) |
| BAM              | BAM file path |

### Options
| workflow | type      | long               | Details                                                                                                                                                                                           |
|----------|-----------|--------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|          | tsv file  | --input            | Path to a tab separated file of sample information and data file paths.                                                                                                                           |
| PON      | bed file  | --pon_bed          | Path to target bed for CNVkit (either capture regions for exome or bins for genome)                                                                                                               |
|          | dir path  | --pon_dir          | Path to directory containing CNVkit panel of normals                                                                                                                                              |
| PON      | string    | --pon_sample       | Name of the sample to be used as repesentative to generate the CNVKit target bed file.                                                                                                            |
| PON      | tsv file  | --pon_tsv          | Path to a tab separated file of normal coverage file paths.                                                                                                                                       |


### Folder locations
By default, Docker containers cannot access files located on the machine they run on. Therefore, local folders need to be mapped to folders inside the container using
``-v local_folder:container_folder`` when calling `docker run`:
```bash
sudo docker run \
--user $(id -u):$(id -g) \
-v <your_working_directory>:/var/pipeline/ \
-v <your_ref_directory>:/var/pipeline/ref/ \
-v <your_temp_directory>:/var/pipeline/temp/ \
rolandradlab/mocaseq:<mocaseq_version> <options>
```
Importantly, the pipeline requires that: \
The main working directory needs to be mapped to ``/var/pipeline/``. \
The reference directory (``ref``, NOT ``ref/GRCm38.p6``) needs to be mapped to ``/var/pipeline/ref/ ``. \
The temp directory needs to be mapped to ``/var/pipeline/temp/ ``.

### Interactive mode
By default, this Docker image automatically runs the MoCaSeq pipeline, as detailed above. However, all scripts included in the ``repository`` folder can be used separately, to allow adjustment of specific aspects of the pipeline. For this, the image can be started in interactive mode:
```bash
sudo docker run \
-it --entrypoint=/bin/bash \
--user $(id -u):$(id -g) \
-v <your_working_directory>:/var/pipeline/ \
rolandradlab/mocaseq:<mocaseq_version>
```

## TL;DR
1. Download and install [Docker](https://www.docker.com/products/docker-desktop).

2. Set the name and location of the working directory. This will be used for testing and the reference files will be located here. Use a volume with at least 250 GB of free disk space.
```bash
working_directory=/PATH/TO/WORKING_DIRECTORY
```

3. Create the working directory:
```bash
mkdir -p ${working_directory} \
&& cd ${working_directory}
```

4. Download and unzip the repository from Github:
```bash
wget https://github.com/roland-rad-lab/MoCaSeq/archive/master.zip \
&& unzip master.zip \
&& rm master.zip \
&& mv MoCaSeq-master ${working_directory}/MoCaSeq
```

5. Download the Docker image from Dockerhub:
```bash
sudo docker pull rolandradlab/mocaseq
```

6. Test the pipeline, which automatically downloads the required reference files. If succesful, reference files will be located in `ref/` and test results in `MoCaSeq_Test/`. This will take up to 24 hours!
```bash
sudo docker run \
-v ${working_directory}:/var/pipeline/ \
rolandradlab/mocaseq \
--test yes
```

7. Use the provided script to download both tumor and matched normal FASTQ files from one pancreatic cancer, which developed in a conditionally-activated Kras<sup>G12D</sup>-model. `all` will download both WES (100x) and WGS (30x) data, using 100 GB of disk space. Use `WES` or `WGS` to only download the respective files.
```bash
mkdir -p ${working_directory}/raw \
&& cd ${working_directory}/raw \
&& sh ${working_directory}/MoCaSeq/repository/Preparation_GetExemplaryData.sh WES \
&& cd ${working_directory}
```
* The raw data is available from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) using the run accessions ERR2230828 (WES Tumour), ERR2230866 (WES Normal), ERR2210078 (WGS Tumour) and ERR2210079 (WGS Normal).

8. Now run test the pipeline using a "real-life" sample. Replace `<threads>` and `<ram>`, then run the pipeline using the data downloaded in Step 7. Depending on the available CPU and RAM, this will take approximately 24 hours.
```bash
sudo docker run \
--user $(id -u):$(id -g) \
-v ${working_directory}:/var/pipeline/ \
rolandradlab/mocaseq \
-nf '/var/pipeline/raw/S821-WES.Normal.R1.fastq.gz' \
-nr '/var/pipeline/raw/S821-WES.Normal.R2.fastq.gz' \
-tf '/var/pipeline/raw/S821-WES.Tumor.R1.fastq.gz' \
-tr '/var/pipeline/raw/S821-WES.Tumor.R2.fastq.gz' \
--name S821-WES \
--sequencing_type WES \
--threads <threads> \
--RAM <RAM> \
--artefact GT
```

9. Browse the results located in `S821/results/`.

## Bug reports
Please send comments and bug reports to: sebastian.lange [@] tum.de

## Citation
Please cite the corresponding protocol published concurrently to this repository: \
\
S. Lange, T. Engleitner, S. Mueller, R. Maresch, M. Zwiebel, L. Gonzaléz-Silva, G. Schneider, R. Banerjee, F. Yang, G. Vassiliou, M. Friedrich, D. Saur, I. Varela and R. Rad (2020).
Analysis pipelines for cancer genome sequencing in mice. *Nature Protocols*. \
https://doi.org/10.1038/s41596-019-0234-7

Specific versions of the repository can be cited directly: \
\
S. Lange. MoCaSeq: Analysis pipelines for cancer genome sequencing in mice. \
https://doi.org/10.5281/zenodo.3344535

## License
The program is distributed under the MIT license.
