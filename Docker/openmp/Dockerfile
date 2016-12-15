FROM buzmakov/baseimage
MAINTAINER Alexey Buzmakov <buzmakov@gmail.com>

ADD wpg_install.sh /opt/wpg_install.sh

RUN DEBIAN_FRONTEND=noninteractive  apt-get install -y --no-install-recommends  \
		wget \
		build-essential \
		unzip \
		libgomp1 &&\

 bash /opt/wpg_install.sh &&\
 apt-get remove -y wget build-essential unzip &&\
 apt-get clean && \
 apt-get autoremove -y && apt-get autoclean -y &&\
 rm -rf /var/lib/apt/lists/*

VOLUME /data
WORKDIR /data
EXPOSE 8888
ENV PYTHONPATH=$PYTHONPATH:/opt/WPG-feature-openmp
RUN cd /opt/WPG-feature-openmp/ && python -m pytest . && md5sum -c tests/reference.md5

CMD jupyter notebook --no-browser --ip 0.0.0.0 --notebook-dir /