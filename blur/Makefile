# Author: David Holmqvist <daae19@student.bth.se>

CXX       = g++
CXXFLAGS  = -std=c++17 -g -Wunused -Wall -Wunused
# If you hit pthread link issues on your platform, try -lpthread instead of -pthread:
PTHREADS  = -pthread

.PHONY: all clean

all: blur blur_par

# ----------------- Sequential (original) -----------------

blur: matrix ppm filters blur.cpp
	$(CXX) $(CXXFLAGS) blur.cpp matrix.o ppm.o filters.o -o blur

filters: matrix filters.hpp filters.cpp
	$(CXX) $(CXXFLAGS) -c filters.cpp -o filters.o

matrix: matrix.hpp matrix.cpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

ppm: ppm.hpp ppm.cpp
	$(CXX) $(CXXFLAGS) -c ppm.cpp -o ppm.o

# ----------------- Parallel (pthreads) -----------------

# Sources you added from my parallelization step:
#   - filters_par.hpp
#   - filters_par.cpp
#   - blur_par.cpp

filters_par.o: filters_par.cpp filters_par.hpp matrix.hpp
	$(CXX) $(CXXFLAGS) -c filters_par.cpp -o filters_par.o

blur_par: blur_par.cpp matrix.o ppm.o filters_par.o
	$(CXX) $(CXXFLAGS) blur_par.cpp matrix.o ppm.o filters_par.o -o blur_par $(PTHREADS)

# ----------------- Housekeeping -----------------

clean:
	rm -rf blur blur_par *.ppm *.o *.dSYM 2> /dev/null
