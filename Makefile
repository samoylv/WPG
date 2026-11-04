# WPG make targets
#
#   make install                  official SRW wheel (srwpy from PyPI) + WPG in editable mode
#   make install SRW_MODE=omp     SRW built from sources with OpenMP + WPG (Linux)
#   make srw SRW_MODE=cuda        only SRW, built from sources with CUDA (Linux)
#
# SRW_MODE: wheel (default) | cpu | omp | cuda, see build/build_srw.py --help
# SRW_REF:  SRW branch, tag or commit for source builds (default: master)
# SRW_SRC:  existing SRW checkout for source builds (instead of cloning to build/SRW)
# SRW_CUDA_ARCHS: native (GPU of this machine, default) | all (80 86 89 90 120)
# PYTHON:   Python of the target environment (default: python3)

PYTHON ?= python3
SRW_MODE ?= wheel
SRW_REF ?= master
SRW_SRC ?=
SRW_CUDA_ARCHS ?= native
SRW_ARGS = --mode $(SRW_MODE) --ref $(SRW_REF) --python $(shell $(PYTHON) -c "import sys; print(sys.executable)") --cuda-archs $(SRW_CUDA_ARCHS) $(if $(SRW_SRC),--src $(SRW_SRC))

all: install

srw:
ifeq ($(SRW_MODE),wheel)
	$(PYTHON) -m pip install --upgrade "srwpy>=4.2"
else
	$(PYTHON) build/build_srw.py $(SRW_ARGS)
endif

install: srw
	$(PYTHON) -m pip install -e ".[test]"

test:
	$(PYTHON) -m pytest

test-notebooks:
	$(PYTHON) -m pytest -m notebooks

lint:
	$(PYTHON) -m ruff check .

doc:
	cd docs; make html

clean:
	rm -rf build/SRW build/srw_build_info.json

.PHONY: all srw install test test-notebooks lint doc clean
