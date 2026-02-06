CXXFLAGS = -std=c++11 -pthread -std=c++17 -m64 -I/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.28.08/x86_64-ubuntu22-gcc114-opt/include 
LDFLAGS = -L/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.28.08/x86_64-ubuntu22-gcc114-opt/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.28.08/x86_64-ubuntu22-gcc114-opt/lib -pthread -lm -ldl -rdynamic
EIGEN_PATH = /usr/include/eigen3 
