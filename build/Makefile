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
	wget --no-check-certificate -nc -O $(SRC_DIR)/srw.zip https://github.com/buzmakov/SRW/archive/feature/srw_lite_openmp.zip; echo 'Done'
	unzip -o $(SRC_DIR)/srw.zip -d $(TMP_DIR)
	rm -rf $(TMP_DIR)/SRW
	mv $(TMP_DIR)/SRW-feature-srw_lite_openmp $(TMP_DIR)/SRW
	patch $(TMP_DIR)/SRW/cpp/py/Makefile py3.patch
	cd $(TMP_DIR)/SRW; make all;
	mkdir -p $(BUILD_DIR)/lib
	cp $(TMP_DIR)/SRW/cpp/gcc/srwlpy*.so $(BUILD_DIR)/lib/
	cp $(TMP_DIR)/SRW/cpp/gcc/srwlpy*.so $(TMP_DIR)/SRW/env/work/srw_python
