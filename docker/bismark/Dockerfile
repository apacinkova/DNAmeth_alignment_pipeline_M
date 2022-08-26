# Bismark: bisulfite sequencing alignment

# base image: Ubuntu
FROM ubuntu:20.04 AS build-base

# software versions
ENV seqtk_version=1.3
ENV samtools_version=1.15.1
ENV bowtie2_version=2.4.5
ENV bismark_version=0.23.0
ENV DEBIAN_FRONTEND=noninteractive


FROM build-base AS build-seqtk

WORKDIR /
RUN apt-get update --fix-missing \
    && apt-get -y install make gcc zlib1g-dev

# install seqtk
ADD https://github.com/lh3/seqtk/archive/v${seqtk_version}.tar.gz v${seqtk_version}.tar.gz
RUN tar -xzvf v${seqtk_version}.tar.gz \
    && rm v${seqtk_version}.tar.gz \
    && mv seqtk-${seqtk_version} seqtk \
    && cd seqtk \
    && make


FROM build-base AS build-samtools

WORKDIR /
RUN apt-get update --fix-missing \
    && apt-get -y install make gcc zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev

# install samtools
ADD https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 ${samtools_version}.tar.bz2
RUN tar xjvf ${samtools_version}.tar.bz2 \
    && rm ${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --prefix=/samtools \
    && make \
    && make install


FROM build-base AS build-bowtie

WORKDIR /
RUN apt-get update --fix-missing \
    && apt-get -y install unzip libncurses5-dev

# install bowtie2
ADD https://github.com/BenLangmead/bowtie2/releases/download/v${bowtie2_version}/bowtie2-${bowtie2_version}-linux-x86_64.zip bowtie2-${bowtie2_version}-linux-x86_64.zip
RUN unzip bowtie2-${bowtie2_version}-linux-x86_64.zip \
    && rm bowtie2-${bowtie2_version}-linux-x86_64.zip \
    && mv bowtie2-${bowtie2_version}-linux-x86_64 bowtie2


FROM build-base AS build-bismark

WORKDIR /
RUN apt-get update --fix-missing \
    && apt-get -y install libgd-graph-perl

# install bismark
ADD https://github.com/FelixKrueger/Bismark/archive/refs/tags/${bismark_version}.tar.gz .
RUN tar -xzvf ${bismark_version}.tar.gz \
    && mv Bismark-${bismark_version} Bismark


FROM ubuntu:20.04
RUN apt-get update && apt-get install -y \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \    
    perl
COPY --from=build-seqtk /seqtk/* /usr/local/bin/
COPY --from=build-samtools /samtools/bin/* /usr/local/bin/
COPY --from=build-bowtie /bowtie2/bowtie2 /bowtie2/bowtie2-* /usr/local/bin/
COPY --from=build-bismark /Bismark/*bismark* /Bismark/coverage2cytosine /Bismark/bam2nuc /Bismark/filter_non_conversion /Bismark/methylation_consistency /usr/local/bin/
