Introduction
========
This repository serves as a companion to an analysis worklow protocol for mouse cancer next-generation data.

Abstract
========
Mouse models of human cancer have proven invaluable in linking genetics to mechanisms and phenotypes. The tremendous opportunities for reverse as well as forward genetics in mice gain further momentum through the sequencing revolution. Sequencing analysis pipelines were however developed for humans, and hence do not account for species-specific differences in genome structures and experimental setups.

Here, we describe standardised computational pipelines tailored specifically for the analysis of mouse genomic data. We present workflows for the detection of all genetic alterations, including single nucleotide variants, indels, copy number variation, loss of heterozygosity and complex rearrangements. All components have been extensively validated and cross-compared using multiple methodologies. The protocol also contributes novel analytical tools, such as pipelines for inference of chromothripsis.

We provide code for all workflows, give step by step guidance on the conduction of individual analysis types and provide advice for data interpretation. The protocol takes 2-7 days, depending on the desired analyses.

Installation
========
1. This requires a number of tools, used for the detection of somatic mutations. Installation procedures are listed in `repository/Preparation_SystemSetup.sh`.

2. Edit `config.sh` to update the relevant paths for each tool, as well as the temporary directory and the number of CPU threads and RAM to be used.

3. Run `sh repository/Preparation_GetReferenceData.sh config.sh` to automatically download all reference files needed to the current directory (inside a newly created folder `Genomes`).

Usage
========
The complete pipeline is wrapped inside a shell-script. During the analysis, the .fastq-files are automatically copied to the target directory.
```
sh MouseCancerGenomeAnalysis.sh <name> \
<fastq_normal_1> <fastq_normal_2> <fastq_tumor_1> <fastq_tumor_2> \
<sequencing_type> <config_file> [none|CT|GT] [phred33|phred64]"
```

Example
========

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


Citation
========
Please cite our primary paper: \
Mueller, S., Engleitner, T., Maresch, R., Zukowska, M., Lange, S., …, Rad, R. (2018). \
Evolutionary routes and KRAS dosage define pancreatic cancer phenotypes. Nature, 554(7690), 62–68. \
https://doi.org/10.1038/nature25459