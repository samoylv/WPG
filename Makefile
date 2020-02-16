all: core modules
OPENMP_MODE ?= 0

core:
	cd build; OPENMP_MODE=$(OPENMP_MODE) make srw

modules:
	cp build/lib/srwlpy*.so wpg/srw/
	cp build/tmp/SRW/cpp/src/lib/srwlib.h docs/

clean:
	cd build; make clean
	rm wpg/srw/srwlpy*.so
	rm docs/srwlib.h

doc:
	cd docs; make html

test:
	pytest

.PHONY: all core modules clean doc test
