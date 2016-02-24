all: core modules
	
core:
	cd build; make srw

modules: core
	cp build/lib/srwlpy*.so wpg/
	cp build/tmp/SRW/env/work/srw_python/*.py wpg/
	cp -r build/tmp/SRW/env/work/srw_python/data_example_* wpg/
	cp build/tmp/SRW/cpp/src/lib/srwlib.h docs/

clean:
	cd build; make clean
	rm wpg/srwlpy*.so
	rm wpg/srwlib.py
	rm wpg/uti_plot.py
	rm docs/srwlib.h

ipython: core
	cd modules/ipython_env; make all
	cd modules; make ipython
	cp lib/srwlpy*.so samples/srw_python
	cp lib/srwlpy*.so samples/srw_python/wavefront
	cd samples/srw_python; ln -s ../../modules/ipython_env/build/bin/ipython

doc:
	cd docs; make html
	
.PHONY: all core modules clean ipython doc
