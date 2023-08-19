FROM ubuntu:18.04

LABEL org.label-schema.name="MoCaSeq"
LABEL org.label-schema.vcs-url="https://github.com/roland-rad-lab/MoCaSeq"

LABEL author="niklas.andrade-kraetzig [@] tum.de"
LABEL maintainer="niklas.andrade-kraetzig [@] tum.de"
LABEL org.label-schema.url="https://www.imo.med.tum.de"

ENV TARGET_DIR /var/pipeline
ENV PACKAGE_DIR=/opt
ENV TEMP_DIR=/var/tmp
ENV HOME /var

## adds the folder with local installation files to the global envir
ADD /data/resources /data/resources

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

ENV LANGUAGE=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

ENV PATH ${PATH}:${PACKAGE_DIR}/bin

ENV PERL5LIB ${PACKAGE_DIR}/vep-96:${PERL5LIB}

## Configure default locale, mostly to avoid problems with R (sorting etc.),
## see https://github.com/rocker-org/rocker/issues/19
## define a user to be available via runtime flag '--user docker'
## install essential system packages
RUN	apt update \
	&& apt install -y --no-install-recommends \
		locales \
	&& sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
	&& dpkg-reconfigure --frontend=noninteractive locales \
	&& update-locale LANG=en_US.UTF-8 \
	&& mkdir -p ${PACKAGE_DIR}/bin \
	&& apt upgrade -y \
	&& apt install -y --no-install-recommends \
		apt-utils \
		ca-certificates \
		dirmngr \
		fonts-texgyre \
		gpg-agent \
		software-properties-common \
		gosu \ 
                libharfbuzz-dev \
                libfribidi-dev \
	&& apt -y --no-install-recommends upgrade \
	&& apt install -y --no-install-recommends \
		apt-transport-https \
		alien \
		ant \
		autoconf \
		automake \
		bioperl \
		build-essential \
		bowtie \
		bc \
		cmake \
		cpanminus \
		curl \
		default-jre \
		g++ \
		gdebi \
		git \
		git-all \
		git-daemon-sysvinit \
		gparted \
		gnupg2 \
		libbam-dev \
		libbz2-dev \
		libcairo2-dev \
		libcurl3-dev \
		libexpat1-dev \
		libgd-dev \
		libgsl0-dev \
		libjemalloc-dev \
		liblapacke-dev \
		liblzma-dev \
		libncurses5-dev \
		libpng-dev \
		libssl-dev \
		libxml2-dev \
		openssh-server \
		pandoc \
		parallel \
		perl \
		perl-base \
		poppler-utils \
		python \
		python-dev \
		rpm \
		rsync \
		sra-toolkit \
		unzip \
		wget \
		xvfb \
		zip \
		zlib1g-dev \
	&& add-apt-repository 'ppa:openjdk-r/ppa' -y \
	&& add-apt-repository 'ppa:git-core/ppa' -y \
	&& apt update \
	&& apt install -y --no-install-recommends \
		openjdk-8-jdk \
	&& rm -rf /var/lib/apt/lists/* \
	&& cd ${TEMP_DIR} \
	&& curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash \
	&& apt-get install git-lfs \
	&& git lfs install

RUN 	curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py \
	&& python2 get-pip.py \
	&& rm get-pip.py \
	&& pip2 install numpy \
	&& pip2 install scipy \
	matplotlib \
	six \
	deeptools \
	pandas \
	multiprocessing \
	PyVCF \
	ConfigParser \
	Cheetah \
	pysam \
	fisher

RUN	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 \
	&& apt-get update \
	&& add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
	&& apt update \
	&& apt install -y --no-install-recommends \
		r-base=3.6.3-1bionic \
		r-base-dev=3.6.3-1bionic \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get update \
	&& apt-get install -y --no-install-recommends libgit2-dev

RUN apt-get -y install libpoppler-cpp-dev

RUN R -e 'install.packages("usethis")'
RUN R -e 'install.packages("devtools")'
RUN R -e 'devtools::install_github("mskcc/facets", build_vignettes = F)'
RUN R -e 'devtools::install_github("mskcc/pctGCdata")'
RUN R -e 'install.packages("numDeriv")'
RUN R -e 'install.packages("/data/resources/ABSOLUTE_1.0.6.tar.gz",repos=NULL,type="source")'
RUN R -e 'install.packages(pkgs=c("BiocManager"),dependencies=TRUE)'
RUN R -e 'install.packages("XML", repos = "http://www.omegahat.net/R")'
RUN R -e 'install.packages("bit")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz",repos=NULL)'
RUN R -e 'devtools::install_version("devEMF", version = "4.0", repos = "http://cran.us.r-project.org")'
RUN R -e 'install.packages("rapportools")'
RUN R -e 'install.packages("pdftools")'
RUN R -e 'install.packages("ggpubr")'
RUN R -e 'install.packages("plotly")'
RUN R -e 'install.packages("htmlwidgets")'

RUN R -e 'BiocManager::install(pkgs=c("tidyverse","splitstackshape","GenomicRanges","optparse","zoo","ggplot2","CopywriteR","HMMcopy","DNAcopy","GenomeInfoDb","Biostrings","data.table","RColorBrewer","pheatmap","biomaRt","BSgenome.Mmusculus.UCSC.mm10","BSgenome.Hsapiens.UCSC.hg38","BSgenome.Hsapiens.UCSC.hg19","deconstructSigs","LSD","randtests","svglite","dupRadar","TitanCNA","doMC","naturalsort","SomaticSignatures","SomaticCancerAlterations", "BubbleTree", "SNPchip"),version="3.10",ask=FALSE,update=TRUE)'

RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz",repos=NULL)' # package ff is updated in biocmanager, but we need this version

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
	&& apt-get update \
	&& apt install -y --no-install-recommends python3.7 python3-pip

RUN apt-get install -y --no-install-recommends python-numpy python-scipy python-matplotlib python-reportlab python-pandas
RUN apt-get install -y --no-install-recommends python3.7-dev

RUN apt install -y --no-install-recommends python3.7 \
	&& curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
	&& python3.7 get-pip.py \
	&& rm get-pip.py

RUN pip3.7 install multiqc
RUN pip3.7 install pomegranate==0.14.5 # only this version works for cnvkit
RUN pip3.7 install cnvkit==0.9.9

RUN apt-get install -y --no-install-recommends python-tk \
	&& python3.7 -m pip install seaborn \
	&& python3.7 -m pip install sinfo \
	&& python3.7 -m pip install datatable \
	&& python3.7 -m pip install scanpy \
	&& python3.7 -m pip install pysam \
	&& python3.7 -m pip install dask \
	&& python3.7 -m pip install toolz \
	&& python3.7 -m pip install 'fsspec>=0.3.3' \
	&& python3.7 -m pip install openpyxl

# SAMBLASTER
RUN cd ${TEMP_DIR} \
	&& git clone https://github.com/GregoryFaust/samblaster.git \
	&& cd samblaster \
	&& git checkout 'v.0.1.26' \
	&& make \
	&& cd .. \
	&& mv samblaster ${PACKAGE_DIR}/samblaster-0.1.26/

RUN cd ${TEMP_DIR} \
	&& curl -fsSL https://download.docker.com/linux/debian/gpg | apt-key add - \
	&& add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable edge" \
	; apt update \
	&& apt install -y --no-install-recommends \
	docker-ce docker-ce-cli containerd.io

RUN	cd ${TEMP_DIR} \
# bwa-mem v0.7.17 (http://bio-bwa.sourceforge.net)
	&& cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2' \
	&& tar -xjvf bwa-0.7.17.tar.bz2 \
	&& cd bwa-0.7.17 \
	&& make \
	&& cd ${TEMP_DIR} \
	&& mv ./bwa-0.7.17 ${PACKAGE_DIR}/ \
	&& ln -sf ${PACKAGE_DIR}/bwa-0.7.17/bwa ${PACKAGE_DIR}/bin/bwa \
	&& rm -rf bwa-0.7.17.tar.bz2 \
# temporary files from bwa.kit are needed for human ref genome generation
	&& cd ${TEMP_DIR} \
	&& wget -nv 'https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2' \
	&& tar -xjf bwakit-0.7.15_x64-linux.tar.bz2 \
	&& mv ./bwa.kit/resource-GRCh38 ${PACKAGE_DIR}/bwa-0.7.17/bwakit/ \
	&& mv ./bwa.kit/resource-human-HLA ${PACKAGE_DIR}/bwa-0.7.17/bwakit/ \
	&& mv bwa.kit ${PACKAGE_DIR}/bwakit-0.7.15 \
	&& rm bwakit-0.7.15_x64-linux.tar.bz2

# htslib 1.10.2 and bcftools 1.10.2 and samtools 1.10
RUN	cd ${TEMP_DIR} \
	&& git clone https://github.com/samtools/htslib.git \
	&& git clone https://github.com/samtools/samtools.git \
	&& git clone https://github.com/samtools/bcftools.git \
	&& cd htslib \
	&& git checkout '1.10.2' \
	&& make \
	&& make install \
	&& cd ../bcftools \
	&& git checkout '1.10.2' \
	&& make \
	&& make install \
	&& cd ../samtools \
	&& git checkout '1.10' \
	&& make \
	&& make install \
	&& cd ${TEMP_DIR} \
	&& rm -rf htslib bcftools samtools

# bedtools v2.28.0 (https://github.com/arq5x/bedtools2)
RUN	cd ${TEMP_DIR} \
	&& git clone https://github.com/arq5x/bedtools2.git \
	&& cd bedtools2 \
	&& git checkout 'v2.28.0' \
	&& make \
	&& make install \
	&& cd ${TEMP_DIR} \
	&& rm -rf bedtools2

# Delly2 v0.8.3 (https://github.com/dellytools/delly)
RUN	apt update \
	&& apt install -y --no-install-recommends \
	libboost-date-time-dev \
	libboost-program-options-dev \
	libboost-system-dev \
	libboost-filesystem-dev \
	libboost-iostreams-dev \
	&& cd ${TEMP_DIR} \
	&& wget https://github.com/dellytools/delly/releases/download/v0.8.3/delly_v0.8.3_linux_x86_64bit \
	&& mkdir ${PACKAGE_DIR}/delly-0.8.3/ \
	&& mv delly_v0.8.3_linux_x86_64bit ${PACKAGE_DIR}/delly-0.8.3/delly \
	&& chmod a+x ${PACKAGE_DIR}/delly-0.8.3/delly \
	&& ln -sf ${PACKAGE_DIR}/delly-0.8.3/delly ${PACKAGE_DIR}/bin/delly

# Fasta-to-Fastq (https://github.com/ekg/fasta-to-fastq/)
RUN	cd ${TEMP_DIR} \
	&& git clone https://github.com/ekg/fasta-to-fastq.git \
	&& cp ./fasta-to-fastq/fasta_to_fastq.pl ${PACKAGE_DIR}/bin/ \
	&& rm -rf fasta-to-fastq

# GATK v3.8.1.0 (https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive\&version=3.8-1-0-gf15c1c3ef)
RUN	cd ${TEMP_DIR} \
	&& wget -nv -O gatk-3.8.1.0.tar.bz2 'https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2' \
	&& tar -xvjf gatk-3.8.1.0.tar.bz2 \
	&& mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/ ${PACKAGE_DIR}/gatk-3.8.1.0 \
	&& rm gatk-3.8.1.0.tar.bz2

# GATK v4.1.7.0 (https://software.broadinstitute.org/gatk/)
RUN	cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/broadinstitute/gatk/releases/download/4.1.7.0/gatk-4.1.7.0.zip' \
	&& unzip gatk-4.1.7.0.zip \
	&& mkdir -p ${PACKAGE_DIR}/gatk-4.1.7.0 \
	&& cp ./gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ${PACKAGE_DIR}/gatk-4.1.7.0/ \
	&& ln -sf ${PACKAGE_DIR}/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ${PACKAGE_DIR}/gatk-4.1.7.0/gatk.jar \
	&& rm -rf gatk-4.1.7.0.zip gatk-4.1.7.0 \
# GATK v4.1.4.1 (https://software.broadinstitute.org/gatk/)
	&& cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip' \
	&& unzip gatk-4.1.4.1.zip \
	&& mkdir -p ${PACKAGE_DIR}/gatk-4.1.4.1 \
	&& cp ./gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar ${PACKAGE_DIR}/gatk-4.1.4.1/ \
	&& ln -sf ${PACKAGE_DIR}/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar ${PACKAGE_DIR}/gatk-4.1.4.1/gatk.jar \
	&& rm -rf gatk-4.1.4.1.zip gatk-4.1.4.1 \
# GATK v4.1.3.0 (https://software.broadinstitute.org/gatk/)
	&& cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/broadinstitute/gatk/releases/download/4.1.3.0/gatk-4.1.3.0.zip' \
	&& unzip gatk-4.1.3.0.zip \
	&& mkdir -p ${PACKAGE_DIR}/gatk-4.1.3.0 \
	&& cp ./gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar ${PACKAGE_DIR}/gatk-4.1.3.0/ \
	&& ln -sf ${PACKAGE_DIR}/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar ${PACKAGE_DIR}/gatk-4.1.3.0/gatk.jar \
	&& rm -rf gatk-4.1.3.0.zip gatk-4.1.3.0 \
# GATK v4.1.0.0 (https://software.broadinstitute.org/gatk/)
	&& cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/broadinstitute/gatk/releases/download/4.1.0.0/gatk-4.1.0.0.zip' \
	&& unzip gatk-4.1.0.0.zip \
	&& mkdir -p ${PACKAGE_DIR}/gatk-4.1.0.0 \
	&& cp ./gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar ${PACKAGE_DIR}/gatk-4.1.0.0/ \
	&& ln -sf ${PACKAGE_DIR}/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar ${PACKAGE_DIR}/gatk-4.1.0.0/gatk.jar \
	&& rm -rf gatk-4.1.0.0.zip gatk-4.1.0.0 \
# GATK v4.2.0.0 (https://software.broadinstitute.org/gatk/)
	&& cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip' \
	&& unzip gatk-4.2.0.0.zip \
	&& mkdir -p ${PACKAGE_DIR}/gatk-4.2.0.0 \
	&& cp ./gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar ${PACKAGE_DIR}/gatk-4.2.0.0/ \
	&& ln -sf ${PACKAGE_DIR}/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar ${PACKAGE_DIR}/gatk-4.2.0.0/gatk.jar \
	&& rm -rf gatk-4.2.0.0.zip gatk-4.2.0.0 \
# GATK requires Java 8 to be active
	&& update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java

# HMMCopy Utils (https://github.com/shahcompbio/hmmcopy_utils/)
RUN	cd ${TEMP_DIR} \
	&& git clone https://github.com/shahcompbio/hmmcopy_utils.git \
	&& cd hmmcopy_utils \
	&& cmake . \
	&& make \
	&& mkdir -p ${PACKAGE_DIR}/hmmcopy_utils \
	&& cp -R ./bin ./util ${PACKAGE_DIR}/hmmcopy_utils/ \
	&& cd ${TEMP_DIR} \
	&& rm -rf hmmcopy_utils \
# fastQC v0.11.8 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
	&& cd ${TEMP_DIR} \
	&& wget -nv 'www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip' \
	&& unzip fastqc_v0.11.8.zip \
	&& mv FastQC ${PACKAGE_DIR}/FastQC-0.11.8 \
	&& chmod +x ${PACKAGE_DIR}/FastQC-0.11.8/fastqc \
	&& ln -sf ${PACKAGE_DIR}/FastQC-0.11.8/fastqc ${PACKAGE_DIR}/bin/fastqc \
	&& rm -rf fastqc_v0.11.8.zip
# Picard v2.20.0 (https://broadinstitute.github.io/picard)
RUN	cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/broadinstitute/picard/releases/download/2.20.0/picard.jar' \
	&& mkdir -p ${PACKAGE_DIR}/picard-2.20.0 \
	&& mv picard.jar ${PACKAGE_DIR}/picard-2.20.0

# sambamba v0.7.0 (https://broadinstitute.github.io/picard)
RUN	cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/biod/sambamba/releases/download/v0.7.0/sambamba-0.7.0-linux-static.gz' \
	&& gunzip sambamba-0.7.0-linux-static.gz \
	&& mkdir -p ${PACKAGE_DIR}/sambamba-0.7.0 \
	&& mv sambamba-0.7.0-linux-static ${PACKAGE_DIR}/sambamba-0.7.0/sambamba-0.7.0 \
	&& chmod +x ${PACKAGE_DIR}/sambamba-0.7.0/sambamba-0.7.0 \
	&& ln -sf ${PACKAGE_DIR}/sambamba-0.7.0/sambamba-0.7.0 ${PACKAGE_DIR}/bin/sambamba

# SnpEff v4.3T (http://snpeff.sourceforge.net)
RUN	cd ${TEMP_DIR} \
	&& wget -nv 'https://netix.dl.sourceforge.net/project/snpeff/snpEff_v4_3t_core.zip' \
	&& unzip 'snpEff_v4_3t_core.zip' \
	&& mv clinEff ${PACKAGE_DIR}/clinEff-4.3T \
	&& mv snpEff ${PACKAGE_DIR}/snpEff-4.3T \
	&& rm snpEff_v4_3t_core.zip \
	&& wget -nv -O snpEff_v4_3_GRCh38.92.zip 'https://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh38.92.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsnpeff%2Ffiles%2Fdatabases%2Fv4_3%2FsnpEff_v4_3_GRCh38.92.zip%2Fdownload&ts=1556478741' \
	&& unzip snpEff_v4_3_GRCh38.92.zip -d snpEff_v4_3_GRCh38.92 \
	&& mkdir -p ${PACKAGE_DIR}/snpEff-4.3T/data \
	&& mv snpEff_v4_3_GRCh38.92/data/GRCh38.92 ${PACKAGE_DIR}/snpEff-4.3T/data/GRCh38.92 \
	&& echo 'GRCh38.92.genome : Homo_sapiens' >> ${PACKAGE_DIR}/snpEff-4.3T/snpEff.config \
	&& echo 'GRCh38.92.reference : ftp://ftp.ensembl.org/pub/release-92/gtf/' >> ${PACKAGE_DIR}/snpEff-4.3T/snpEff.config \
	&& rm snpEff_v4_3_GRCh38.92.zip \
	&& rm -r snpEff_v4_3_GRCh38.92 \
	&& wget -nv -O snpEff_v4_3_GRCm38.86.zip 'https://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCm38.86.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsnpeff%2Ffiles%2Fdatabases%2Fv4_3%2FsnpEff_v4_3_GRCm38.86.zip%2Fdownload&ts=1557205135' \
	&& unzip snpEff_v4_3_GRCm38.86.zip -d snpEff_v4_3_GRCm38.86 \
	&& mv snpEff_v4_3_GRCm38.86/data/GRCm38.86 ${PACKAGE_DIR}/snpEff-4.3T/data/GRCm38.86 \
#	&& echo 'GRCm38.86.genome : Homo_sapiens' >> ${PACKAGE_DIR}/snpEff-4.3T/snpEff.config \
#	&& echo 'GRCm38.86.reference : ftp://ftp.ensembl.org/pub/release-86/gtf/' >> ${PACKAGE_DIR}/snpEff-4.3T/snpEff.config \
	&& rm snpEff_v4_3_GRCm38.86.zip \
	&& rm -r snpEff_v4_3_GRCm38.86

# Trimmomatic v0.39 (http://www.usadellab.org)
RUN	cd ${TEMP_DIR} \
	&& wget -nv 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip' \
	&& unzip Trimmomatic-0.39.zip \
	&& mv Trimmomatic-0.39 ${PACKAGE_DIR}/trimmomatic-0.39 \
	&& rm Trimmomatic-0.39.zip

# msisensor v0.5 (https://github.com/ding-lab/msisensor, MSI testing)
RUN	cd ${TEMP_DIR} \
	&& git clone 'https://github.com/ding-lab/msisensor.git' \
	&& cd msisensor \
	&& git checkout '0.5' \
	&& make \
	&& cp msisensor ${PACKAGE_DIR}/bin \
	&& cd ${TEMP_DIR} \
	&& rm -rf msisensor

# Bam-matcher
RUN	cd ${TEMP_DIR} \
	&& git clone 'https://bitbucket.org/sacgf/bam-matcher.git' \
	&& mv ./bam-matcher ${PACKAGE_DIR}/

# Strelka 2.9.10 (https://github.com/Illumina/strelka)
RUN	cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2' \
	&& tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2 \
	&& mv ./strelka-2.9.10.centos6_x86_64 ${PACKAGE_DIR}/strelka-2.9.10 \
	&& rm strelka-2.9.10.centos6_x86_64.tar.bz2

# Manta 1.6.0 (https://github.com/Illumina/manta)
RUN	cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2' \
	&& tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2 \
	&& mv manta-1.6.0.centos6_x86_64 ${PACKAGE_DIR}/manta-1.6.0 \
	&& rm manta-1.6.0.centos6_x86_64.tar.bz2

# VariantQC 1.07 (https://github.com/BimberLab/DISCVRSeq)
RUN	cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/BimberLab/DISCVRSeq/releases/download/1.07/DISCVRSeq-1.07.jar' \
	&& mkdir -p ${PACKAGE_DIR}/DISCVRSeq-1.07/ \
	&& mv DISCVRSeq-1.07.jar ${PACKAGE_DIR}/DISCVRSeq-1.07/

# Vcftools 0.1.16 (https://github.com/vcftools/vcftools.git)
RUN	cd ${TEMP_DIR} \
	&& git clone 'https://github.com/vcftools/vcftools.git' \
	&& cd ./vcftools \
	&& git checkout 'v0.1.16' \
	&& ./autogen.sh \
	&& ./configure \
	&& make \
	&& make install \
	&& cd ${TEMP_DIR} \
	&& rm -rf vcftools

# new cmake version to fix bamtools installation
RUN cd ${TEMP_DIR} \
	&& mkdir cmake \
	&& wget https://github.com/Kitware/CMake/releases/download/v3.21.0-rc2/cmake-3.21.0-rc2.tar.gz \
	&& tar -zxvf cmake-3.21.0-rc2.tar.gz \
	&& cd cmake-3.21.0-rc2 \
	&& ./bootstrap \
	&& make \
	&& make install \
	&& hash -r

RUN apt-get install libjsoncpp-dev

# BAMTOOLS GetBaseCounts v.1.2.2 (https://github.com/zengzheng123/GetBaseCountsMultiSample.git)
RUN	cd ${TEMP_DIR} \
	&& git clone https://github.com/pezmaster31/bamtools.git \
	&& cd bamtools/ \
	&& mkdir -p build \
	&& cd build/ \
	&& cmake -DCMAKE_INSTALL_PREFIX=${TEMP_DIR}/bamtools .. \
	&& make \
	&& make install \
	&& cd ${TEMP_DIR} \
	&& git clone https://github.com/zengzheng123/GetBaseCountsMultiSample.git \
	&& cd GetBaseCountsMultiSample \
	&& git checkout 'v1.2.2' \
	&& g++ -o3 -I${TEMP_DIR}/bamtools/include/bamtools -L${TEMP_DIR}/bamtools/lib/ GetBaseCountsMultiSample.cpp -lbamtools -lz -o GetBaseCountsMultiSample -fopenmp \
	&& cp ./GetBaseCountsMultiSample ${PACKAGE_DIR}/bin/ \
	&& cd ${TEMP_DIR} \
	&& rm -rf bamtools GetBaseCountsMultiSample

# Ensembl VEP 96.0 (https://github.com/Ensembl/ensembl-vep.git)
RUN	cd ${TEMP_DIR} \
	&& git clone 'https://github.com/Ensembl/ensembl-vep.git' \
	&& cd ensembl-vep \
	&& git checkout 'release/96.0' \
	&& cd ${TEMP_DIR} \
	&& mv ensembl-vep ${PACKAGE_DIR}/vep-96 \
	&& perl ${PACKAGE_DIR}/vep-96/INSTALL.pl --AUTO a --DESTDIR ${PACKAGE_DIR}/vep-96 --NO_UPDATE --NO_TEST --NO_HTSLIB

# vcf2maf 1.6.17 (https://github.com/mskcc/vcf2maf/archive/v1.6.17.tar.gz)
RUN	cpanm --notest LWP::Simple Archive::Zip Archive::Extract HTTP::Tiny Test::Simple File::Copy::Recursive Perl::OSType Module::Metadata version TAP::Harness CGI Encode CPAN::Meta JSON DBD::SQLite Set::IntervalTree Archive::Tar Time::HiRes Module::Build Bio::Root::Version \
	&& cd ${TEMP_DIR} \
	&& wget -nv 'https://github.com/mskcc/vcf2maf/archive/v1.6.17.tar.gz' \
	&& tar -xzf v1.6.17.tar.gz \
	&& mv ./vcf2maf-1.6.17 ${PACKAGE_DIR} \
	&& rm -rf v1.6.17.tar.gz

# this is part of FACETS in R
RUN g++ -std=c++11 /usr/local/lib/R/site-library/facets/extcode/snp-pileup.cpp -lhts -o /usr/local/lib/R/site-library/facets/extcode/snp-pileup

WORKDIR /var/pipeline/

RUN cd ${PACKAGE_DIR} \
	&& git clone --branch human-pipeline https://github.com/roland-rad-lab/MoCaSeq.git \
	&& cd ${PACKAGE_DIR}/MoCaSeq/ \
    && chmod 775 entrypoint.sh \
    && chmod 775 MoCaSeq.sh \
    && chmod 775 repository/Meta_logstats.sh

ENTRYPOINT ["/opt/MoCaSeq/entrypoint.sh"]
CMD ["-h"]
