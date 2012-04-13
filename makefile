# Change the following to match your environment.
# Specifically, you may have to add -Lpath or -Ipath to FLAGS to point to correct paths of boost/blas/lapack libraries.
CXX=g++
CXXFLAGS=-DNDEBUG -O3 -march=native -Wall -c
CC=gcc
CFLAGS=-O0 -Wall -c
FC=gcc
FCFLAGS=-m64 -c
LDFLAGS=-lblas -llapack -lpthread -lgfortran -lrt


CXXSOURCES=src/data.cpp src/data_model.cpp src/main.cpp src/matrix.cpp src/model.cpp src/prior.cpp src/sampler.cpp src/symmmatrix.cpp src/utils.cpp src/vector.cpp src/inih/cpp/INIReader.cpp
CSOURCES=src/inih/ini.c
OBJECTS=$(CXXSOURCES:.cpp=.cpp.o) $(CSOURCES:.c=.c.o) src/dchex.f.o
EXECUTABLE=bmagwa

all: src/*.hpp $(CXXSOURCES) $(CSOURCES) src/dchex.f $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

%.cpp.o: %.cpp src/*.hpp # force recompilation if any headers change
	$(CXX) $(CXXFLAGS) $< -o $@
%.c.o: %.c
	$(CC) $(CFLAGS) $< -o $@
%.f.o: %.f
	$(FC) $(FCFLAGS) $< -o $@

.PHONY: clean
clean:
	rm $(OBJECTS) $(EXECUTABLE)
