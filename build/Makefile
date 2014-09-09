ORIGDIR = $(CURDIR)
BUILD_DIR = $(ORIGDIR)
BIN_DIR = $(BUILD_DIR)/bin
TMP_DIR = $(ORIGDIR)/tmp
SRC_DIR = $(ORIGDIR)/sources
PYTHON_TOMO = $(BIN_DIR)/python
EASY_INSTALL = $(BIN_DIR)/easy_install
PIP = $(BIN_DIR)/pip
CFLAGS += -fPIC
CPPFLAGS += -fPIC

PY_VERSION = 2.7.6

.PHONY: all numpy matplotlib ipython modules pip python setup_tools clean pyyaml srw

all: python pip ipython

python: $(BIN_DIR)/python

$(BIN_DIR)/python: $(SRC_DIR)/Python-$(PY_VERSION).tgz
	mkdir -p $(TMP_DIR)
	#mkdir -p $(BUILD_DIR)
	mkdir -p $(SRC_DIR)
		
	tar xzvf $(SRC_DIR)/Python-$(PY_VERSION).tgz --directory=$(TMP_DIR)
	cd $(TMP_DIR)/Python-$(PY_VERSION); ./configure --prefix=$(BUILD_DIR) --with-zlib; make -j8 all; make install;

$(SRC_DIR)/Python-$(PY_VERSION).tgz:
	mkdir -p $(SRC_DIR)
	cd $(SRC_DIR); wget -nc http://www.python.org/ftp/python/$(PY_VERSION)/Python-$(PY_VERSION).tgz

setup_tools: python $(BIN_DIR)/easy_install

$(BIN_DIR)/easy_install: $(SRC_DIR)/setuptools-0.6c11.tar.gz
	tar xzvf $(SRC_DIR)/setuptools-0.6c11.tar.gz --directory=$(TMP_DIR)
	cd $(TMP_DIR)/setuptools-0.6c11/;$(PYTHON_TOMO) setup.py install;cd $(ORIGDIR)

$(SRC_DIR)/setuptools-0.6c11.tar.gz:
	cd $(SRC_DIR); wget -nc http://pypi.python.org/packages/source/s/setuptools/setuptools-0.6c11.tar.gz


ipython: modules srw
	$(PIP) install pyzmq
	$(PIP) install tornado
	$(PIP) install jinja2
	$(PIP) install ipython
	rm ../gcc/srwlpy.so
	cd tmp/SRW-master/cpp/py; $(BUILD_DIR)/bin/pythonbin/python setup.py build_ext --build-lib='../gcc';
	cp $(TMP_DIR)/SRW-master/cpp/gcc/srwlpy.so $(TMP_DIR)/SRW-master/env/work/srw_python
	cd samples; ln -s ../build/bin/ipython

# pyyaml: python
# 	#cd $(SRC_DIR); wget http://pyyaml.org/download/libyaml/yaml-0.1.4.tar.gz
# 	#tar xzvf $(SRC_DIR)/yaml-0.1.4.tar.gz --directory=$(TMP_DIR)
# 	#cd $(TMP_DIR)/yaml-0.1.4; ./configure --prefix=$(BUILD_DIR); make install;

# 	cd $(SRC_DIR); wget http://pyyaml.org/download/pyyaml/PyYAML-3.10.tar.gz
# 	tar xzvf $(SRC_DIR)/PyYAML-3.10.tar.gz --directory=$(TMP_DIR)
# 	cd $(TMP_DIR)/PyYAML-3.10; $(PYTHON_TOMO) setup.py --without-libyaml install

modules: python
	$(PIP) install numpy
	$(PIP) install matplotlib
	$(PIP) install h5py


pip: setup_tools
	$(EASY_INSTALL) pip

clean:
	rm -rf $(TMP_DIR)
	rm -rf $(BUILD_DIR)/include
	rm -rf $(BUILD_DIR)/info
	rm -rf $(BUILD_DIR)/lib
	rm -rf $(BUILD_DIR)/tmp


srw:
	mkdir -p $(TMP_DIR)
	
	cd $(SRC_DIR); wget -nc http://www.fftw.org/fftw-2.1.5.tar.gz
	tar xzvf $(SRC_DIR)/fftw-2.1.5.tar.gz --directory=$(TMP_DIR)
	cd $(TMP_DIR)/fftw-2.1.5; ./configure  --prefix=$(BUILD_DIR) --with-pic --enable-float; make clean; make -j ;make install;

	wget --no-check-certificate -nc -O $(SRC_DIR)/srw.zip https://github.com/ochubar/SRW/archive/ebdb80779f333edcdb35e3cd966869cbaf9538c0.zip; echo 'Done'
	unzip -o $(SRC_DIR)/srw.zip -d $(TMP_DIR)
	rm -rf $(TMP_DIR)/SRW
	mv -u $(TMP_DIR)/SRW-ebdb80779f333edcdb35e3cd966869cbaf9538c0 $(TMP_DIR)/SRW
	cp $(BUILD_DIR)/lib/libfftw.a $(TMP_DIR)/SRW/ext_lib/
	cp $(BUILD_DIR)/lib/librfftw.a $(TMP_DIR)/SRW/ext_lib/

	cd $(TMP_DIR)/SRW; make all;
	cp $(TMP_DIR)/SRW/cpp/gcc/libsrw.a $(BUILD_DIR)/lib
	cp $(TMP_DIR)/SRW/cpp/gcc/srwlpy.so $(BUILD_DIR)/lib
	cp $(TMP_DIR)/SRW/cpp/gcc/srwlpy.so $(TMP_DIR)/SRW/env/work/srw_python