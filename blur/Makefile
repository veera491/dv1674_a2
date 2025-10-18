# Author: David Holmqvist <daae19@student.bth.se>

CXX       = g++
CXXFLAGS  = -std=c++17 -O3 -g -Wall -Wunused
PTHREADS  = -pthread

.PHONY: all clean

all: blur blur_par

# --- Sequential ---
blur: blur.cpp matrix.o ppm.o filters.o
	$(CXX) $(CXXFLAGS) blur.cpp matrix.o ppm.o filters.o -o blur

filters.o: filters.cpp filters.hpp matrix.hpp
	$(CXX) $(CXXFLAGS) -c filters.cpp -o filters.o

matrix.o: matrix.cpp matrix.hpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

ppm.o: ppm.cpp ppm.hpp matrix.hpp
	$(CXX) $(CXXFLAGS) -c ppm.cpp -o ppm.o

# --- Parallel (pthreads) ---
filters_par.o: filters_par.cpp filters_par.hpp matrix.hpp
	$(CXX) $(CXXFLAGS) -c filters_par.cpp -o filters_par.o

blur_par: blur_par.cpp matrix.o ppm.o filters_par.o
	$(CXX) $(CXXFLAGS) blur_par.cpp matrix.o ppm.o filters_par.o -o blur_par $(PTHREADS)

clean:
	rm -rf blur blur_par *.ppm *.o *.dSYM 2> /dev/null
