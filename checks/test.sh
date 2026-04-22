root << EOF
  TFile *_file0 = TFile::Open("../Original/data/EmptyFileTestBeamPhase2.root")
  gSystem->AddIncludePath("-I/usr/include/eigen3/");
  .L cubic_fit.C++g
  PulseShape->cd()
  resample(grPulseShape)
EOF
