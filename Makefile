# ===============================
# LoCo C++ Project Makefile
# Dependencies: nanoflann & igraph, boost, gfortran, arpack
# ===============================

# Directories
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin
TOOL_DIR = tools
INC_DIR = $(SRC_DIR)
LIB = src/lib
SUBDIRS := $(shell find $(INC_DIR) -type d)
INCLUDE_DIRS := $(addprefix -I,$(SUBDIRS))

#add boost disr to include flags - important for windows and macOS
INCLUDE_DIRS += $(if $(BOOST_INCLUDE),-I$(BOOST_INCLUDE),)

#add path to nanoflann
#IG_INCLUDE += -I dependencies/nanoflann/include
NANO_INCLUDE = -Iinst/include

CXXFLAGS = -std=c++17 -O3 -Wall -Wextra $(INCLUDE_DIRS) $(NANO_INCLUDE)
LDFLAGS  = -lboost_program_options -lz
# add LTO only for Linux/Mac
ifneq ($(IS_LINUX),)
	CXXFLAGS += -flto=5
else ifneq ($(IS_DARWIN),)
	CXXFLAGS += -flto
endif

# Source files
SRC_FILES := \
  $(SRC_DIR)/SCParser.cpp \
  $(SRC_DIR)/GraphData.cpp \
  $(SRC_DIR)/GraphHandler.cpp \
  $(SRC_DIR)/Neighborhood.cpp \
  $(TOOL_DIR)/LoCo.cpp

OBJ_FILES := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))


#######################################
# PLATFORM DETECTION
#######################################
#definitions for systems
UNAME_S := $(shell uname -s)

IS_LINUX  := $(filter Linux,$(UNAME_S))
IS_DARWIN := $(filter Darwin,$(UNAME_S))
IS_WIN    := $(filter MINGW% MSYS% CYGWIN%,$(UNAME_S))

#######################################
# BOOST INSTALL PATHS
#######################################

#######################################
# macOS (Homebrew)
#######################################

ifeq ($(IS_DARWIN),Darwin)
    BOOST_PREFIX := $(shell brew --prefix boost 2>/dev/null)

    BOOST_CPPFLAGS += -I$(BOOST_PREFIX)/include
    BOOST_LDFLAGS  += -L$(BOOST_PREFIX)/lib
endif

#######################################
# Linux
#######################################

ifneq ($(IS_LINUX),)
    BOOST_CPPFLAGS += $(shell pkg-config --cflags boost-program-options 2>/dev/null)
    BOOST_LDFLAGS  += $(shell pkg-config --libs boost-program-options 2>/dev/null)
endif

#######################################
# Windows (vcpkg)
#######################################

ifneq ($(IS_WIN),)
    VCPKG_ROOT ?= C:/vcpkg
    BOOST_TRIPLET := x64-mingw-static

    BOOST_CPPFLAGS += -I$(VCPKG_ROOT)/installed/$(BOOST_TRIPLET)/include
    BOOST_LDFLAGS  += -L$(VCPKG_ROOT)/installed/$(BOOST_TRIPLET)/lib
endif

#######################################
# FINAL BOOST FLAGS (USED BY COMPILER)
#######################################

CXXFLAGS += $(BOOST_CPPFLAGS)
LDFLAGS  += $(BOOST_LDFLAGS)
LDLIBS   += $(BOOST_LIBS)

#######################################
# INSTALL DEPENDENCIES (SYSTEM-DEPENDENT)
#######################################
install:
# -------------------------------
# NANOFLANN is included as a header only in the repository - therefore we origionally used the branch below
# -------------------------------	
#	mkdir -p inst/include
#	cd inst/include && \
#	if [ ! -d nanoflann ]; then \
#		git clone https://github.com/jlblancoc/nanoflann --branch v1.3.2; \
#	else \
#		echo "nanoflann already exists, skipping clone"; \
#	fi

#######################################
# LINUX (APT / SYSTEM)
#######################################
	@if [ "$(IS_LINUX)" = "Linux" ]; then \
		echo "Installing dependencies (Linux)..."; \
		sudo apt-get update && sudo apt-get install -y \
			libboost-program-options-dev \
			zlib1g-dev; \
	fi

#######################################
# macOS (Homebrew)
#######################################
	@if [ "$(IS_DARWIN)" = "Darwin" ]; then \
		echo "Installing dependencies (macOS)..."; \
		brew update; \
		brew install boost zlib || true; \
	fi

#######################################
# Windows (vcpkg)
#######################################
	@if echo "$(UNAME_S)" | grep -E -q "MINGW|MSYS|CYGWIN"; then \
		echo "Installing dependencies (Windows via vcpkg)..."; \
		$(VCPKG_ROOT)/vcpkg install \
			boost-program-options \
			zlib \
			--triplet x64-mingw-static; \
	fi

	mkdir -p $(BUILD_DIR)
	mkdir -p $(BIN_DIR)
	@echo "Dependencies installed."

# ===============================
# Build LoCo
# ===============================

loco: $(OBJ_FILES) | $(BIN_DIR)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJ_FILES) $(LDFLAGS) $(LDLIBS) -o $(BIN_DIR)/loco

# ===============================
# Compile sources to object files
# ===============================
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

#$(BUILD_DIR)/LoCo.o: tools/LoCo.cpp
#	@mkdir -p $(dir $@)
#	$(CXX) $(CXXFLAGS) -c $< -o $@

# ===============================
# Directories
# ===============================
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# ===============================
# Clean
# ===============================
.PHONY: clean install loco
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

#test dataset has 4 correlated var that range through grpah from low to high corr, 4 medium constant corr and a bunch of non corr
#this dataset is a bit smaller and noisy, we do not check results automatically
test_loco_a:
	./bin/loco -i ./test/simulatedData1.tsv -o bin -p test_a -n 20 -x 0.3 -s 25 -t 1 -q 1

#simple test where we have 50 non corr variables of two clsuters, then in each cluster diff 3 var correlated
#you should run it with printing cliques and find in each neighborhood roughly ONLY 1,2,3 or 4,5,6
# RESULT: this should always find the two correlatuion sets A,B,C and D,E,F
test_loco_b:
	./bin/loco -i ./test/simulatedData2.tsv -o bin -p test_b -c -v ./test/cellStateGenes.txt -w ./test/cellSignalGenes.txt -t 10
	./test/test_b.sh

#tets with 1000 cells and 58 features
# 20 totally reandom
# 5 (A,B,C,D,E) of correaltions (between all the same) from 0.1 up to 1 (in 100 steps with 10cells in every 'step')
#additioannly 20 of corr 0.75 (same, no change in corr). However, we order the cells to follow the increasing corr of first 4 AB (A, B,C ,D,E)
# So the ABs of low correlation are also low in 1:20-0.75_corr and as those increase (but same corr) the corr of A,B,C,D,E increases
#=> idea: const corr of 0.75 gives graph the structure and following this structure the corr of first 4 should increase
#THIS SHOULD DETECT ONLY THE FIRST 5 ABs, the other 20 have good correlations, but those DO NOT CHANGE along the cell-manifold
# we run 5X: ./bin/loco -i ./test/simulatedData3.tsv -o bin -p test_c -c -n 20 -s 100 -x 0.5 -t 50 -q 2 -m 2 -a 0.01
# and check that 4times we have in the top 5 correlations any pairs of A,B,C,D,E & this set is reported everytime and we have p-values<0.05
test_loco_c:
	./test/test_c.sh
#increasing n to 100 makes everything significant

# simulate the signlaing markers also as markers with a signoidal activation
test_run_loco_sigmoidal:
	/usr/bin/time ./bin/loco -i ./test/data_1.tsv -o bin -p data_1 -c -n 100 -s 50 -x 0.4 -z 1 -t 50 -m 2 -q 2 -a 0.01 -u 1000 -f 0
# test aboce script 5 times and make sure in the top 5 correlations pairs we only see pairs of the middle program (Ms) or an end program (Es)
test_loco_sigmoidal:
	./test/test_sigmoidal.sh

test_loco_sigmoidal_granularities:
	/usr/bin/time ./bin/loco -i ./test/data_1.tsv -o bin -p data_1 -c -n 100 -s [10,50,100,200] -x 0.4 -z 1 -t 50 -m 2 -q 2 -a 0.01 -u 1000 -f 0

#-v ./test/paperCellstateMarkers.txt
#-v ./test/paperCellstateMarkers.txt -w ./test/paperCellsignalMarkers.txt

# simulate the signaling markers as uniformal distributions
test_loco_uniform:
	/usr/bin/time ./bin/loco -i ./test/data_2.tsv -o bin -p data_2 -c -n 1000 -s 100 -x 0.5 -z 1 -t 10 -f 20 -v test/paperCellstateMarkers.txt -w test/paperCellsignalMarkers.txt

#5K cells,m when having more than 50N p-values seem to not make sense anymore
test_loco_uniform_noSignalMarkers:
	/usr/bin/time ./bin/loco -i ./test/data_2.tsv -o bin -p data_2_b -c -n 50 -s 100 -x 0.5 -z 1 -t 10 -f 20 -v test/paperCellstateMarkers.txt

test_loco:
	make test_loco_a
	make test_loco_b
	make test_loco_c
	make test_loco_sigmoidal
	make test_run_loco_sigmoidal
	make test_loco_uniform
