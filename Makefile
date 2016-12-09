all: core modules
	
core:
	cd build; make srw

modules: core
	cp build/lib/srwlpy*.so wpg/
	cp build/tmp/SRW/cpp/src/lib/srwlib.h docs/

clean:
	cd build; make clean
	rm wpg/srwlpy*.so
	rm docs/srwlib.h

doc:
	cd docs; make html
	
.PHONY: all core modules clean ipython doc
