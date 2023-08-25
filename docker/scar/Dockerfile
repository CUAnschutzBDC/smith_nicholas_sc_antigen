# Get base image
#FROM python:3.9.17-bullseye
FROM continuumio/miniconda3

# environment
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH $PATH:/opt/conda/bin

RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

COPY environment.yaml .
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda install -n base conda-libmamba-solver && \
    conda config --set solver libmamba && \
    conda env update -n base -f environment.yaml && \
    conda clean -afy 
