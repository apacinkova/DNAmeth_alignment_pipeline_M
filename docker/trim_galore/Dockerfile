# TrimGalore! v0.6.6
# https://github.com/FelixKrueger/TrimGalore

# base image: Ubuntu
FROM ubuntu:20.04

# cutadapt version
ARG cutadapt_version=4.1
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /
RUN apt-get update --fix-missing \
    && apt-get -y install --no-install-recommends --no-install-suggests \
        curl perl default-jre python3 unzip build-essential libnss-sss\
        less zlib1g-dev libdeflate-dev libisal-dev\
    && yes '' | apt-get install python3-pip \
    && apt-get install pigz

# install fastqc
ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip /tmp/fastqc.zip
RUN unzip /tmp/fastqc.zip -d /opt/ && \
    rm /tmp/fastqc.zip && \
    chmod u+x /opt/FastQC/fastqc
ENV PATH="/opt/FastQC/:${PATH}"

# install cutadapt
RUN python3 -m pip install --user --upgrade cutadapt==${cutadapt_version} \
    && cp /root/.local/bin/cutadapt /usr/local/bin

# install trim galore
ADD https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz ./0.6.6.tar.gz
RUN tar -xzvf 0.6.6.tar.gz \
    && rm -rf 0.6.6.tar.gz \
    && mv TrimGalore-0.6.6/trim_galore /usr/local/bin

RUN apt-get remove -y unzip \
    && apt-get autoremove -y \
    && apt-get clean

WORKDIR $HOME