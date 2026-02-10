# Makefile using config.mk
include config.mk

CXX = g++ -g -O0

all: CreateData.x multifit.x

CreateData.x: CreateData.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

multifit.x: multifit.cc PulseChiSqSNNLS.cc
	$(CXX) $(CXXFLAGS) -I $(EIGEN_PATH) -o $@ $^ $(LDFLAGS)

clean:
	rm -f CreateData.x multifit.x
