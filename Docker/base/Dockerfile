FROM continuumio/miniconda3
MAINTAINER Alexey Buzmakov <buzmakov@gmail.com>

RUN conda install -y nomkl \
					numpy \
					scipy \
					pyqt=4.11 \
					matplotlib \
					h5py \
					pytest \
					jupyter \
					libgcc &&\
	conda clean -y --all

RUN DEBIAN_FRONTEND=noninteractive  apt-get install -y --no-install-recommends  \
		texlive-fonts-recommended \
        texlive-latex-base \
        texlive-generic-recommended \
        texlive-latex-extra \
        pandoc && \
 apt-get remove -y wget build-essential unzip &&\
 apt-get clean && \
 apt-get autoremove -y && apt-get autoclean -y &&\
 rm -rf /var/lib/apt/lists/*
