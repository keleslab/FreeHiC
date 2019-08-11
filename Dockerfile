FROM ubuntu 
FROM python:3.6
MAINTAINER Ye Zheng<yezheng@stat.wisc.edu>
RUN apt-get update && apt-get install -y samtools bwa bedtools wget git build-essential zlib1g-dev libbz2-dev liblzma-dev
RUN pip install numpy scipy pysam bx-python Cython
RUN git clone --recursive git://github.com/yezhengSTAT/FreeHiC
WORKDIR FreeHiC
