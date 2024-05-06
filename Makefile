all: core modules
OPENMP_MODE ?= omp

core:
	#cd build; OPENMP_MODE=$(OPENMP_MODE) make srw

modules:
	cp /Users/twguest/miniconda3/pkgs/srwpy-4.1.0-py312hede676d_0/lib/python3.12/site-packages/srwpy/srwlpy*.so wpg/srw/


clean:
	cd build; make clean
	rm wpg/srw/srwlpy*.so

doc:
	cd docs; make html

test:
	pytest

.PHONY: all core modules clean doc test
