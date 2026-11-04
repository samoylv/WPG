# WPG with SRW, built from the repository root:
#
#   docker build -f Docker/Dockerfile -t wpg .                          # official SRW wheel
#   docker build -f Docker/Dockerfile --build-arg SRW_MODE=omp -t wpg:omp .
#   docker build -f Docker/Dockerfile --build-arg SRW_MODE=cuda \
#       --build-arg BASE=nvidia/cuda:13.2.1-devel-ubuntu24.04 -t wpg:cuda .
#
# See Docker/readme.md for running notebooks and scripts.

ARG BASE=ubuntu:24.04
FROM ${BASE}

# wheel | cpu | omp | cuda (see build/build_srw.py)
ARG SRW_MODE=wheel
ARG SRW_REF=master

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg \
    PATH=/opt/venv/bin:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-venv python3-dev make git ca-certificates \
        $([ "$SRW_MODE" = "wheel" ] || echo gcc g++ libgomp1) \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/venv

COPY . /opt/WPG
WORKDIR /opt/WPG

# No GPU during "docker build": compile CUDA code for all supported architectures
RUN make install PYTHON=/opt/venv/bin/python SRW_MODE=$SRW_MODE SRW_REF=$SRW_REF SRW_CUDA_ARCHS=all \
    && /opt/venv/bin/python -m pip install --no-cache-dir jupyter \
    && rm -rf build/SRW/env/python/build

VOLUME /data
WORKDIR /data
EXPOSE 8888
CMD ["jupyter", "notebook", "--no-browser", "--ip", "0.0.0.0", "--allow-root", "--notebook-dir", "/"]
