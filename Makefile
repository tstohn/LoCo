
# ===============================
# LoCo C++ Project Makefile
# dependencies: nanoflann & igraph
# ===============================


INC_DIR = src
LIB = src/lib
SUBDIRS := $(shell find $(INC_DIR) -type d)
INCLUDE_DIRS := $(addprefix -I,$(SUBDIRS))

CXXFLAGS = --std=c++17 -c -O3 -Wall -Wextra $(INCLUDE_DIRS) \
	-I $(HOME)/libraries/igraph_libs/include \
	-I dependencies
#possibly add -Wl, rl flags to find lib also during compile time (look at harris)
#TODO: at the moment hard coded paths for libraries needed (as on harris): we need to active scorr_env
LDFLAGS = -L$(HOME)/libraries/igraph_libs/lib \
	-I dependencies \
	-ligraph -larpack -lblas -lgfortran  \
	-lpthread -lz -lboost_program_options -lboost_iostreams

# Source files
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin
SRC_FILES := $(SRC_DIR)/LoCo.cpp \
             $(SRC_DIR)/lib/Parser/SCParser.cpp \
             $(SRC_DIR)/lib/SCGraph/Constructors/GraphData.cpp \
             $(SRC_DIR)/lib/SCGraph/Constructors/GraphHandler.cpp \
             $(SRC_DIR)/lib/SCGraph/Constructors/Neighborhood.cpp

# Object files in build/
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))


# ===============================
# Install dependencies
# ===============================
install:
#get nanoflann
	mkdir -p dependencies
	cd dependencies && \
	if [ ! -d nanoflann ]; then \
		git clone https://github.com/jlblancoc/nanoflann --branch v1.3.2; \
	else \
		echo "nanoflann already exists, skipping clone"; \
	fi
#create bin
	mkdir -p build;
	mkdir -p bin;
#install igraph
	@OS=$(OS); \
	if [ "$$OS" = "Windows_NT" ]; then \
		echo "Windows detected"; \
		echo "Please install igraph manually (e.g. via vcpkg or conda)"; \
	elif [ "$$(uname -s)" = "Darwin" ]; then \
		echo "macOS detected"; \
		brew install igraph; \
	elif [ "$$(uname -s)" = "Linux" ]; then \
		echo "Linux detected"; \
		if command -v apt >/dev/null; then \
			apt install -y libigraph-dev; \
		elif command -v dnf >/dev/null; then \
			dnf install -y igraph-devel; \
		elif command -v pacman >/dev/null; then \
			pacman -S --noconfirm igraph; \
		else \
			echo "Unknown package manager. Install igraph manually."; \
		fi; \
	else \
		echo "Unknown OS. Install igraph manually."; \
	fi


# ===============================
# Build LoCo executable
# ===============================
loco: $(OBJ_FILES) | $(BIN_DIR)
	$(CXX) $(OBJ_FILES) -o $(BIN_DIR)/loco $(LDFLAGS)

# ===============================
# Compile source files to object files
# ===============================
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# ===============================
# Ensure build and bin directories exist
# ===============================
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# ===============================
# Clean build artifacts
# ===============================
.PHONY: clean install loco all
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

#test dataset has 4 correlated var that range through grpah from low to high corr, 4 medium constant corr and a bunch of non corr
test_loco_a:
	./bin/loco -i ./test/simulatedData1.tsv -o bin -p test_a -n 20 -x 0.3 -s 50 -t 1

#simple test where we have 50 non corr variables of two clsuters, then in each cluster diff 3 var correlated
#you should run it with printing cliques and find in each neighborhood roughly ONLY 1,2,3 or 4,5,6
test_loco_b:
	./bin/loco -i ./test/simulatedData2.tsv -o bin -p test_b -c -v ./test/cellStateGenes.txt -w ./test/cellSignalGenes.txt -t 10

#tets with 1000 cells and 58 features
# 20 totally reandom
# 5 (A,B,C,D,E) of correaltions (between all the same) from 0.1 up to 1 (in 100 steps with 10cells in every 'step')
#additioannly 20 of corr 0.75 (same, no change in corr). However, we order the cells to follow the increasing corr of first 4 AB (A, B,C ,D,E)
# So the ABs of low correlation are also low in 1:20-0.75_corr and as those increase (but same corr) the corr of A,B,C,D,E increases
#=> idea: const corr of 0.75 gives graph the structure and following this structure the corr of first 4 should increase
#THIS SHOULD DETECT ONLY THE FIRST 5 ABs, the other 20 have good correlations, but those DO NOT CHANGE along the cell-manifold
test_loco_c:
	./bin/loco -i ./test/simulatedData3.tsv -o bin -p test_c -c -n 10 -s 100 -x 0.5 -t 50
#increasing n to 100 makes everything significant

# simulate the signlaing markers also as markers with a signoidal activation
test_loco_sigmoidal:
	/usr/bin/time ./bin/loco -i ./test/data_1.tsv -o bin -p data_1 -c -n 100 -s 50 -x 0.4 -z -t 50 -m 2 -q 2 -a 0.01 -u 1000 -f 0

test_loco_sigmoidal_granularities:
	/usr/bin/time ./bin/loco -i ./test/data_1.tsv -o bin -p data_1 -c -n 100 -s [10,50,100,200] -x 0.4 -z -t 50 -m 2 -q 2 -a 0.01 -u 1000 -f 0

#-v ./test/paperCellstateMarkers.txt
#-v ./test/paperCellstateMarkers.txt -w ./test/paperCellsignalMarkers.txt

# simulate the signaling markers as uniformal distributions
test_loco_uniform:
	/usr/bin/time ./bin/loco -i ./test/data_2.tsv -o bin -p data_2 -c -n 1000 -s 50 -x 0.5 -z -t 10 -f 20 -v test/paperCellstateMarkers.txt -w test/paperCellsignalMarkers.txt

#5K cells,m when having more than 50N p-values seem to not make sense anymore
test_loco_uniform_noSignalMarkers:
	/usr/bin/time ./bin/loco -i ./test/data_2.tsv -o bin -p data_2_b -c -n 50 -s 100 -x 0.5 -z -t 10 -f 20 -v test/paperCellstateMarkers.txt

test_loco:
	make test_loco_a
	make test_loco_b
	make test_loco_c
	make test_loco_sigmoidal
	make test_loco_uniform
