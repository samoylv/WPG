# WPG in Docker

Build from the repository root (the build context must be the whole repository):

```bash
docker build -f Docker/Dockerfile -t wpg .                              # official SRW wheel (serial)
docker build -f Docker/Dockerfile --build-arg SRW_MODE=omp -t wpg:omp . # SRW with OpenMP
docker build -f Docker/Dockerfile --build-arg SRW_MODE=cuda \
    --build-arg BASE=nvidia/cuda:13.2.1-devel-ubuntu24.04 -t wpg:cuda . # SRW with CUDA
```

`SRW_REF` selects the SRW branch, tag or commit for source builds (`master` by default).

Run the tests:

```bash
docker run --rm wpg python -m pytest /opt/WPG/tests
docker run --rm --gpus all wpg:cuda python -m pytest /opt/WPG/tests
```

Run jupyter notebook, `~/tmp/data` is a local folder for results:

```bash
docker run -it --rm -p 8088:8888 -v ~/tmp/data:/data wpg
```

and open http://localhost:8088. WPG samples are in `/opt/WPG/samples`.

Run a script from the local folder:

```bash
docker run -it --rm -u `id -u`:`id -g` -v ~/tmp/data:/data wpg python /data/my_script.py
```
