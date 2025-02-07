all: core modules
OPENMP_MODE ?= omp

core:
	cd build; OPENMP_MODE=$(OPENMP_MODE) make srw

modules:
	find build -name "srwlpy*.so" -exec install '{}' wpg/srw/ \;


clean:
	cd build; make clean
	rm wpg/srw/srwlpy*.so

doc:
	cd docs; make html

test:
	pytest

.PHONY: all core modules clean doc test
