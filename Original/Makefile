# Makefile using config.mk

include config.mk

CXX = g++

# Default flags
CXXFLAGS_DEFAULT = $(CXXFLAGS)
CXXFLAGS_DEBUG   = -g -O0 -Wall -Wextra $(CXXFLAGS)

all: CreateData multifit

# Normal build
CreateData: CreateData.cc
	$(CXX) $(CXXFLAGS_DEFAULT) -o $@ $^ $(LDFLAGS)

multifit: multifit.cc PulseChiSqSNNLS.cc
	$(CXX) $(CXXFLAGS_DEFAULT) -I $(EIGEN_PATH) -o $@ $^ $(LDFLAGS)

# Debug build
debug_CreateData: CreateData.cc
	$(CXX) $(CXXFLAGS_DEBUG) -o $@ $^ $(LDFLAGS)

debug_multifit: multifit.cc PulseChiSqSNNLS.cc
	$(CXX) $(CXXFLAGS_DEBUG) -I $(EIGEN_PATH) -o $@ $^ $(LDFLAGS)

debug: debug_CreateData debug_multifit

clean:
	rm -f CreateData multifit debug_CreateData debug_multifit
