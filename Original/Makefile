# Makefile using config.mk
include config.mk

CXX = g++ -g -O0

all: CreateData multifit

CreateData: CreateData.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

multifit: multifit.cc PulseChiSqSNNLS.cc
	$(CXX) $(CXXFLAGS) -I $(EIGEN_PATH) -o $@ $^ $(LDFLAGS)

clean:
	rm -f CreateData multifit
