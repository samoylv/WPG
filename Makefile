all: core modules

core:
	cd build; make srw

modules:
	cp build/lib/srwlpy*.so wpg/
	cp build/tmp/SRW/cpp/src/lib/srwlib.h docs/

clean:
	cd build; make clean
	rm wpg/srwlpy*.so
	rm docs/srwlib.h

doc:
	cd docs; make html

test:
	pytest

.PHONY: all core modules clean doc pytest
