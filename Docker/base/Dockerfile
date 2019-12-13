FROM continuumio/miniconda3:4.7.12

MAINTAINER Alexey Buzmakov <buzmakov@gmail.com>

RUN DEBIAN_FRONTEND=noninteractive  apt-get update && apt-get install -y --no-install-recommends  \
		texlive-fonts-recommended \
        texlive-latex-base \
        texlive-generic-recommended \
        texlive-latex-extra \
        pandoc && \
 apt-get clean && \
 apt-get autoremove -y && apt-get autoclean -y &&\
 rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

ENV PYTHONDONTWRITEBYTECODE=true

RUN conda install -y --file requirements.txt \
	&& conda clean -afy

 ENTRYPOINT ["tini", "-g", "--"]
