FROM ubuntu:14.04

MAINTAINER Alexey Buzmakov <buzmakov@gmail.com>

###################
# Setup environment
###################
ENV ROOT_DIR /data/S2E
RUN mkdir -p $ROOT_DIR


###################
# Install packages
###################
ENV PYTHON_MAJOR 2
ENV PYTHON_MINOR 7
RUN apt-get update && apt-get install -y \
    build-essential \
    python${PYTHON_MAJOR}.${PYTHON_MINOR}-dev \
    unzip python-numpy \
    python-scipy python-h5py \ 
    wget \
    nano

###################
# Install WPG
###################
RUN mkdir -p $ROOT_DIR/packages && cd $ROOT_DIR/packages && \
    wget http://github.com/samoylv/WPG/archive/develop.zip -O wpg-develop.zip && \
    unzip wpg-develop.zip && mv WPG-develop WPG && rm wpg-develop.zip && \
    cd WPG && make all && rm -rf build
    
###################
# Install prop
###################

RUN mkdir -p $ROOT_DIR/modules && cd $ROOT_DIR/modules && \
    wget https://github.com/samoylv/prop/archive/develop.zip -O prop-develop.zip && \
    unzip prop-develop.zip && mv prop-develop prop && rm prop-develop.zip


###################
# Setup directories
###################
RUN cd $ROOT_DIR && \
    ln -s /simS2E/workflow workflow && \
    ln -s /simS2E/config config && \
    ln -s /simS2E/data data && \
    ln -s /simS2E/modules modules && \
    ln -s /simS2E/tmp tmp
