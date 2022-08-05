# Bismark: bisulfite sequencing aligment

# base image: Ubuntu
FROM ubuntu:20.04

# software versions
ARG samtools_version=1.15.1
ARG bowtie2_version=2.4.5
ARG bismark_version=0.23.0
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /
RUN apt-get update --fix-missing

#miniconda3: get miniconda3 version 4.7.12
ADD https://repo.continuum.io/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh miniconda.sh
RUN bash miniconda.sh -b -p /opt/miniconda \
    && rm miniconda.sh
    
ENV PATH="/opt/miniconda/bin:${PATH}"

RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge

# install samtools
RUN conda install -y samtools=${samtools_version}

# install bowtie2
RUN conda install -y bowtie2=${bowtie2_version}

# install bismark
RUN conda install -y bismark=${bismark_version} \
    && conda clean -y --all \
    && rm -f /opt/miniconda.sh

RUN apt-get autoremove -y \
    && apt-get clean

WORKDIR $HOME