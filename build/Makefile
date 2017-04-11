ORIGDIR = $(CURDIR)
BUILD_DIR = $(ORIGDIR)
BIN_DIR = $(BUILD_DIR)/bin
TMP_DIR = $(ORIGDIR)/tmp
SRC_DIR = $(ORIGDIR)/sources
PYTHON_TOMO = $(BIN_DIR)/python

.PHONY: all srw

all: srw

clean:
	rm -rf $(TMP_DIR)
	rm -rf $(BUILD_DIR)/include
	rm -rf $(BUILD_DIR)/info
	rm -rf $(BUILD_DIR)/lib
	rm -rf $(BUILD_DIR)/tmp
	rm -rf $(BUILD_DIR)/sources


srw:
	mkdir -p $(TMP_DIR)
	mkdir -p $(BUILD_DIR)/sources
	#wget --no-check-certificate -nc -O $(SRC_DIR)/srw.zip https://github.com/buzmakov/SRW/archive/feature/srw_lite_openmp.zip; echo 'Done'
	### CFG fix to WPG issue #113
	wget --no-check-certificate -nc -O $(SRC_DIR)/srw.zip https://github.com/buzmakov/SRW/archive/80c89e0961965c9d77d6f2f623a4c3f2c2302a5d.zip; echo 'Done'
	unzip -o $(SRC_DIR)/srw.zip -d $(TMP_DIR)
	rm -rf $(TMP_DIR)/SRW
	#mv $(TMP_DIR)/SRW-feature-srw_lite_openmp $(TMP_DIR)/SRW

	### CFG fix to WPG issue #113
	mv $(TMP_DIR)/SRW-80c89e0961965c9d77d6f2f623a4c3f2c2302a5d $(TMP_DIR)/SRW
	cd $(TMP_DIR)/SRW/ext_lib; tar xf fftw-2.1.5.tar.gz; cd fftw-2.1.5; ./configure --enable-float --with-pic; make -j8; cp fftw/.libs/libfftw.a ../.
	### CFG fix to WPG issue #113

	cd $(TMP_DIR)/SRW; make all;
	mkdir -p $(BUILD_DIR)/lib
	cp $(TMP_DIR)/SRW/cpp/gcc/srwlpy*.so $(BUILD_DIR)/lib/
	cp $(TMP_DIR)/SRW/cpp/gcc/srwlpy*.so $(TMP_DIR)/SRW/env/work/srw_python
