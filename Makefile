all: core modules
	
core:
	cd build; make srw

modules: core
	cp build/lib/srwlpy.so wpg/
	cp build/tmp/SRW-master/env/work/srw_python/srwlib.py wpg/
	cp build/tmp/SRW-master/env/work/srw_python/uti_plot.py wpg/

clean:
	cd build; make clean
	rm wpg/srwlpy.so
	rm wpg/srwlib.py
	rm wpg/uti_plot.py

ipython: core
	cd modules/ipython_env; make all
	cd modules; make ipython
	cp lib/srwlpy.so samples/srw_python
	cp lib/srwlpy.so samples/srw_python/wavefront
	cd samples/srw_python; ln -s ../../modules/ipython_env/build/bin/ipython


.PHONY: all core modules clean ipython
