# TrimGalore!
# https://github.com/FelixKrueger/TrimGalore

# base image: Ubuntu
FROM ubuntu:20.04 AS build-base

# software versions
ENV trimgalore_version=0.6.6
ENV DEBIAN_FRONTEND=noninteractive


FROM build-base AS build-fastqc

WORKDIR /
RUN apt-get update --fix-missing \
    && apt-get -y install unzip libncurses5-dev

# install fastqc
ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip /tmp/fastqc.zip
RUN unzip /tmp/fastqc.zip -d /opt/ && \
    rm /tmp/fastqc.zip && \
    chmod u+x /opt/FastQC/fastqc


FROM build-base AS build-cutadapt

WORKDIR /
RUN apt-get update --fix-missing \
    && apt-get -y install python3 python3-pip


FROM build-base AS build-trimgalore

# install trim galore
ADD https://github.com/FelixKrueger/TrimGalore/archive/${trimgalore_version}.tar.gz ./${trimgalore_version}.tar.gz
RUN tar -xzvf ${trimgalore_version}.tar.gz \
    && rm -rf ${trimgalore_version}.tar.gz \
    && mv TrimGalore-${trimgalore_version}/trim_galore trim_galore


FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV cutadapt_version=2.8-2build1

RUN apt-get update && apt-get install -y \
    cutadapt=${cutadapt_version} \
    python3 perl default-jre 

COPY --from=build-fastqc /opt/FastQC/fastqc /usr/local/bin/
COPY --from=build-trimgalore /trim_galore /usr/local/bin/