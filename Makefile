IBM_FLAGS = -qsmp=omp
GNU_STANDARD = c++11

bgp: COMP = bgxlc++_r
bgp: CFLAGS = $(IBM_FLAGS)
gnu: COMP = g++
gnu: CFLAGS = -fopenmp -std=$(GNU_STANDARD)
polus: COMP = xlc++_r
polus: CFLAGS = $(IBM_FLAGS)

.PHONY: bgp polus gnu clean
bgp: bgp_bin
gnu: gnu_bin
polus: polus_bin

BINARY_NAME = solver
SOURCE_DIR = source
BUILD_DIR = build
HEADER_FILES := $(wildcard $(SOURCE_DIR)/*.h)
OBJECT_FILES := $(BUILD_DIR)/main.o $(BUILD_DIR)/specialops.o

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(HEADER_FILES)
	$(COMP) $(CFLAGS) -c -o $@ $<

bgp_bin: $(OBJECT_FILES)
	$(COMP) $(CFLAGS) -o $(BINARY_NAME) $^

gnu_bin: $(OBJECT_FILES)
	$(COMP) $(CFLAGS) -o $(BINARY_NAME) $^

polus_bin: $(OBJECT_FILES)
	$(COMP) $(CFLAGS) -o $(BINARY_NAME) $^

clean:
	rm -f $(BUILD_DIR)/* solver
