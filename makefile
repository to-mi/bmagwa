# Change the following to match your environment.
# Specifically, you may have to add -Lpath or -Ipath to FLAGS to point to correct paths of boost/blas/lapack libraries.
CXX=g++
CXXFLAGS=-std=c++14 -DNDEBUG -O3 -march=native -Wall -c -I/opt/homebrew/include -I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/gcc/include
CC=gcc
CFLAGS=-O0 -Wall -c
FC=gfortran
FCFLAGS=-m64 -c
LDFLAGS=-L/opt/homebrew/opt/gcc/lib/gcc/13 -L/opt/homebrew/opt/openblas/lib -lopenblas -lpthread -lgfortran


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
