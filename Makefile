CUDA_INCLUDE = -I$(CUDA_ROOT)/include
CUDA_CFLAGS = -arch=sm_30 -use_fast_math
CUDA_LFLAGS = -L$(CUDA_ROOT)/lib64 -lcudart
INCLUDE = -Isrc -I$(RATROOT)/include -I$(RATROOT)/src/stlplus
CFLAGS = -g -O3 --std=c++0x -Wall -Werror -ffast-math -fdiagnostics-show-option $(shell root-config --cflags) $(INCLUDE)
LFLAGS = -L$(RATROOT)/lib -lRATEvent_$(shell uname)-g++ -lMinuit $(shell root-config --libs) -ljsoncpp

CC = g++
NVCC = $(CUDA_ROOT)/bin/nvcc

OBJ_DIR = build
SOURCES = src/util.cpp src/fit.cpp
OBJECTS = $(SOURCES:src/%.cpp=$(OBJ_DIR)/%.o)
EXE = fit

all: pdf counting ll.o $(OBJECTS) $(EXE)

clean:
	-$(RM) build/*.o bin/pdf bin/fit

pdf: $(OBJECTS)
	test -d bin || mkdir bin
	$(CC) -o bin/$@ src/$@.cpp $(OBJ_DIR)/util.o $(CFLAGS) $(LFLAGS)

counting: $(OBJECTS)
	test -d bin || mkdir bin
	$(CC) -o bin/$@ src/$@.cpp $(OBJ_DIR)/util.o $(CFLAGS) $(LFLAGS)

$(OBJ_DIR)/%.o: src/%.cpp
	test -d build || mkdir build
	$(CC) -c -o $@ $< $(CFLAGS) $(LFLAGS)

ifndef CUDA_ROOT
$(error CUDA_ROOT is not set)
endif

$(EXE): $(OBJECTS) ll.o
	test -d bin || mkdir bin
	$(CC) -o bin/$@ $(OBJ_DIR)/$@.o $(OBJ_DIR)/util.o $(OBJ_DIR)/ll.o $(CFLAGS) $(LFLAGS) $(CUDA_LFLAGS)

ll.o:
	test -d build || mkdir build
	$(NVCC) -c -o build/ll.o src/ll.cu $(CUDA_CFLAGS) $(CUDA_LFLAGS)

