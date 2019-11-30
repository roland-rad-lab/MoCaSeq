[![DOI](https://zenodo.org/badge/162331694.svg)](https://zenodo.org/badge/latestdoi/162331694)
[![Docker](https://img.shields.io/badge/Docker-avaible-green.svg)](https://cloud.docker.com/repository/docker/rolandradlab/natureprotocols2018/general)
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
	9 Chair of Translational Cancer Research and Institute for Experimental Cancer Therapy, School of Medicine, Technische Universität München, 81675 Munich, Germany<br>
</sub>

### Contents

* [Overview](#Overview)
    + [Abstract](#abstract)
* [System Requirements](#system-requirements)
* [Input formats](#input-formats)
* [Usage](#usage)
    - [User ID](#user-ID)
    - [Folder locations](#folder-locations)
    - [Interactive mode](#interactive-mode)
    - [Options](#options)
* [TL;DR](#tl;dr)
* [Bug reports](#bug-reports)
* [Citation](#citation)
* [License](./LICENSE)

## Overview
This repository serves as a companion to an analysis workflow protocol for mouse cancer next-generation data.

### Abstract
Mouse models of human cancer have transformed our ability to link genetics, molecular mechanisms and phenotypes. Both reverse and forward genetics in mice are currently gaining momentum through advances in next generation sequencing. Methodologies to analyse sequencing data were however developed for humans, and hence do not account for species-specific differences in genome structures and experimental setups. 

Here, we describe standardised computational pipelines tailored specifically for the analysis of mouse genomic data. We present novel tools and workflows for the detection of different alteration types, including single nucleotide variants, indels, copy number variation, loss of heterozygosity and complex rearrangements, such as chromothripsis. Workflows have been extensively validated and cross-compared using multiple methodologies. We also give step by step guidance on the execution of individual analysis types and provide advice on data interpretation. 

The complete code is available [online](https://github.com/roland-rad-lab/MoCaSeq) and as a [dockerized version](https://cloud.docker.com/u/rolandradlab/repository/docker/rolandradlab/mocaseq).

This protocol takes 2-7 days, depending on the desired analyses.  

## System Requirements
We **strongly** recommend to use the Docker version of this pipeline!

- Using the *Docker* version: Platform of your choice.
- Using the *bash* version: Linux, we run this pipeline under Ubuntu 18.04 LTS.
- Hardware: 
	- Minimum: 8-core processor, 32 GB RAM
	- Optimal (running multiple samples in parallel): 48-core processor, 256 GB RAM, Solid-State Drive
- Disk space: 
	- Reference files: 15 GB
	- Results: 30 GB (WES 100x), 300 GB (WGS)
	- Temporary files during analysis: ~170 GB (WES), ~1000 GB (WGS)

All *bash*, *R* and *python* scripts are directly invokable from `repository`.

## Input formats
The standard input format are FASTQ files produced from modern Illumina sequencers. These should be split into forward and reverse reads for both the tumour and the matched normal sample. BAM files can be processed as well, giving the user the choice of starting after alignment (if mapped to GRCm38) or re-mapping all the raw data (see below).

## Usage
We provide a Docker image containing the complete analysis pipeline in order to simplify deployment and to keep software versions as consistent as possible. The basic commandline to run the dockerized pipeline is as follows:
```
sudo docker run \
-v <your_working_directory>:/var/pipeline/ \
rolandradlab/mocaseq:<mocaseq_version> <options>
```
Options are listed [below](#options). When invoked without options, the container will start, display usage information and shut down.

### User ID
Docker by design runs the container and its contents as user root (UID 1 and GID 1). Persistent directories mounted into the container with the option '-v' therefore are owned by root. If this is undesirable, you can pass the UID and GID of your current user into the container by specifying ``--user $(id -u):$(id -g) `` when calling `docker run`:
```
sudo docker run \
-v <your_working_directory>:/var/pipeline/ \
--user $(id -u):$(id -g) \
rolandradlab/mocaseq:<mocaseq_version> <options>
```
This will run the pipeline with your current UID and GID and set the permissions of the output files accordingly.

### Folder locations
By default, Docker containers cannot access files located on the machine they run on. Therefore, local folders need to be mapped to folders inside the container using 
``-v local_folder:container_folder`` when calling `docker run`:
```
sudo docker run \
--user $(id -u):$(id -g) \
-v <your_working_directory>:/var/pipeline/ \
-v <your_ref_directory>:/var/pipeline/ref/ \
-v <your_temp_directory>:/var/pipeline/temp/ \
rolandradlab/mocaseq:<mocaseq_version> <options>
```
Importantly, the pipeline requires that: \
The main working directory needs to be mapped to ``/var/pipeline/``. \
The reference directory (``ref``, NOT ``ref/GRCm38.p6``) needs to be mapped to ``/var/pipeline/ref/ ``.
The temp directory needs to be mapped to ``/var/pipeline/temp/ ``.

### Interactive mode
By default, this Docker image automatically runs the MoCaSeq pipeline, as detailed above. However, all scripts included in the ``repository`` folder can be used separately, to allow adjustment of specific aspects of the pipeline. For this, the image can be started in interactive mode:
```
sudo docker run \
-it --entrypoint=/bin/bash \
--user $(id -u):$(id -g) \
-v <your_working_directory>:/var/pipeline/ \
rolandradlab/mocaseq:<mocaseq_version>
```

### Options
| short    | long                    | Details                                                                                                                                                                                           |
|-------|--------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -n    | --name             | Name of the sample.                                                                                                                                                                               |
| -nf   | --fastq_normal_fw  | Path to first normal Fastq. Do NOT use if running single-sample tumor only.                                                                                                                       |
| -nr   | --fastq_normal_rev | Path to second normal Fastq. Do NOT use if running single-sample tumor only.                                                                                                                      |
| -tf   | --fastq_tumor_fw   | Path to first tumor fastq. Do NOT use if running single-sample normal only.                                                                                                                       |
| -tr   | --fastq_tumor_rev  | Path to second tumor fastq. Do NOT use if running single-sample normal only.                                                                                                                      |
| -nb   | --bam_normal       | Path to normal BAM. Do NOT use in combination with -nf or -nr. When used, -rb MUST be specified.                                                                                                  |
| -tb   | --bam_tumor        | Path to tumor BAM. Do NOT use in combination with -tf or -tr. When used, --repeat_mapping MUST be set to ‘yes’ or ‘no’.                                                                           |
| -rm   | --repeat_mapping   | If -nb or -tb are specified, determines whether mapping is re-done ('yes') or whether the complete mapping procedure is skipped ('no').                                                           |
| -st   | --sequencing_type  | Set to 'WES' or 'WGS'.                                                                                                                                                                            |
| -c    | --config           | Path to configuration file. Optional.                                                                                                                                                             |
| -qc   | --quality_control  | Determines whether QC is done ('yes') or skipped ('no'). Optional.                                                                                                                                |
| -t    | --threads          | Number of CPU threads. Optional. Defaults to 8.                                                                                                                                                   |
| -r    | --RAM              | Amount of Gb RAM. Optional. Defaults to 32.                                                                                                                                                       |
| -temp | --temp_dir         | Path to temporary directory. Optional. Defaults to current working directory.                                                                                                                     |
| -art  | --artefact         | Set to 'GT' (oxidation artefact), 'CT' (FFPE artefact) or 'none'. Optional. If set to something other than 'none' AND Mutect2 is 'yes', forces quality_control to 'yes’. Defaults to none.        |
| -filt | --filtering        | Set to 'all' (AF >= 0.05, , Variant in Tumor >= 2, Variant in Normal <= 1, Coverage >= 5), 'hard' (AF >= 0.1, , Variant in Tumor >= 3, Variant in Normal = 0, Coverage >= 10) or 'none' (no filters). Optional. Defaults to 'hard'.                 |
| -p    | --phred            | If not set, script will try to automatically extract phred-score. Otherwise, set manually to 'phred33' or 'phred64'. 'phred64' only relevant for Illumina data originating before 2011. Optional. |
| -mu   | --Mutect2          | Set to 'yes' or 'no'. Needed for LOH analysis and Titan. Greatly increases runtime for WGS. Optional. Defaults to 'yes'.                                                                          |
| -de   | --Delly            | Set to 'yes' or 'no'. Needed for chromothripsis inference. Do not use for WES. Optional. Defaults to 'no'. Only use in matched-sample mode.                                                       |
| -ti   | --Titan            | Set to 'yes' or 'no'. Greatly increases runtime for WGS. If set to 'yes’, forces Mutect2 to 'yes'. Optional. Defaults to 'yes' for WES and 'no' for WGS. Only use in matched-sample mode.         |
|       | --test             | If set to 'yes': Will download reference files (if needed) and start a test run.                                                                                                                  |
|       | --memstats         | If integer > 0 specified, will write timestamped memory usage and cumulative CPU time usage of the docker container to ./results/memstats.txt every <integer> seconds.                            |
|       | --help             | Show this help.                                                                                                                                                                                   |

## TL;DR
1. Download and install [Docker](https://www.docker.com/products/docker-desktop).

2. Set the name and location of the working directory. This will be used for testing and the reference files will be located here. Use a volume with at least 250 GB of free disk space.
```
working_directory=/PATH/TO/WORKING_DIRECTORY
```

3. Create the working directory:
```
mkdir -p ${working_directory} \
&& cd ${working_directory}
```

4. Download and unzip the repository from Github: 
```
wget https://github.com/roland-rad-lab/MoCaSeq/archive/master.zip \
&& unzip master.zip \
&& rm master.zip \
&& mv MoCaSeq-master ${working_directory}/MoCaSeq
```

5. Download the Docker image from Dockerhub:
```
sudo docker pull rolandradlab/mocaseq
```

6. Test the pipeline, which automatically downloads the required reference files. If succesful, reference files will be located in `ref/` and test results in `MoCaSeq_Test/`. This will take up to 24 hours!
```
sudo docker run \
-v ${working_directory}:/var/pipeline/ \
rolandradlab/mocaseq \
--test yes
``` 

7. Use the provided script to download both tumor and matched normal FASTQ files from one pancreatic cancer, which developed in a conditionally-activated Kras<sup>G12D</sup>-model. `all` will download both WES (100x) and WGS (30x) data, using 100 GB of disk space. Use `WES` or `WGS` to only download the respective files.
```
mkdir -p ${working_directory}/raw \
&& cd ${working_directory}/raw \
&& sh ${working_directory}/MoCaSeq/repository/Preparation_GetExemplaryData.sh WES \
&& cd ${working_directory}
```
* The raw data is available from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) using the run accessions ERR2230828 (WES Tumour), ERR2230866 (WES Normal), ERR2210078 (WGS Tumour) and ERR2210079 (WGS Normal).

8. Now run test the pipeline using a "real-life" sample. Replace `<threads>` and `<ram>`, then run the pipeline using the data downloaded in Step 7. Depending on the available CPU and RAM, this will take approximately 24 hours.
```
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
This repository can be cited using: \
S. Lange. MoCaSeq: Analysis pipelines for cancer genome sequencing in mice. \
https://doi.org/10.5281/zenodo.3344535

Primary paper for which this pipeline was developed: \
S. Mueller, T. Engleitner, R. Maresch, M. Zukowska, S. Lange, […], R. Rad (2018). \
Evolutionary routes and KRAS dosage define pancreatic cancer phenotypes. Nature, 554(7690), 62–68. \
https://doi.org/10.1038/nature25459

## License
The program is distributed under the MIT license.