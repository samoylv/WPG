FROM continuumio/miniconda3

MAINTAINER Alexey Buzmakov <buzmakov@gmail.com>

COPY requirements.txt .

ENV PYTHONDONTWRITEBYTECODE=true

RUN conda install -y --file requirements.txt \
	&& conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && find /opt/conda/lib/python*/site-packages/bokeh/server/static -follow -type f -name '*.js' ! -name '*.min.js' -delete


RUN DEBIAN_FRONTEND=noninteractive  apt-get update && apt-get install -y --no-install-recommends  \
		texlive-fonts-recommended \
        texlive-latex-base \
        texlive-generic-recommended \
        texlive-latex-extra \
        pandoc && \
 apt-get clean && \
 apt-get autoremove -y && apt-get autoclean -y &&\
 rm -rf /var/lib/apt/lists/*

 ENTRYPOINT ["tini", "-g", "--"]
