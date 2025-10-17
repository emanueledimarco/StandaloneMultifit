//
// MultiFit amplitude reconstruction
// To run:
// > g++ -o Example06 Example06.cc PulseChiSqSNNLS.cc -std=c++11 `root-config --cflags --glibs`
// > ./Example06
//

#include <iostream>
#include "PulseChiSqSNNLS.h"
#include "Pulse.h"

#include "TTree.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TFile.h"

using namespace std;

Pulse pSh;

FullSampleVector fullpulse(FullSampleVector::Zero());
FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());
SampleMatrix noisecor(SampleMatrix::Zero());
BXVector activeBX;
SampleVector amplitudes(SampleVector::Zero());

TFile *fout;
TH1D *h01;



void initHist()
{
  fout = new TFile("output.root","recreate");
  h01 = new TH1D("h01", "dA", 1000, -5.0, 5.0);
}

void init()
{
  initHist();
  
  // intime sample is [2]
  double pulseShapeTemplate[NSAMPLES+2];
  for(int i=0; i<(NSAMPLES+2); i++){
    
//     iwf/4. - (500 / 2) + 25. 
//     double x = double( IDSTART + NFREQ * (i + 3) - WFLENGTH / 2);
//     double x = double( IDSTART + NFREQ * (i + 3) - WFLENGTH / 2);
    double x = double( IDSTART + NFREQ * (i + 3) - 500 / 2); //----> 500 ns is fixed!  
    //     x = double( IDSTART + NFREQ * i + 3*25. - 500 / 2. );  //----> 500 ns is fixed!  
    
    pulseShapeTemplate[i] = pSh.fShape(x);
    std::cout << " >>  pulseShapeTemplate[" << i << "] " <<  pulseShapeTemplate[i] << " at x = " << x << std::endl;
    
  }
  //  for(int i=0; i<(NSAMPLES+2); i++) pulseShapeTemplate[i] /= pulseShapeTemplate[2];
  for (int i=0; i<(NSAMPLES+2); ++i) fullpulse(i+7) = pulseShapeTemplate[i];
  
  
  for (int i=0; i<NSAMPLES; ++i) {
    for (int j=0; j<NSAMPLES; ++j) {
      int vidx = std::abs(j-i);
      noisecor(i,j) = pSh.corr(vidx);
    }
  }
  
  int activeBXs[] = { -5, -4, -3, -2, -1,  0,  1,  2,  3,  4 };
  activeBX.resize(10);
  for (unsigned int ibx=0; ibx<10; ++ibx) {
    activeBX.coeffRef(ibx) = activeBXs[ibx];
  } 
  //  activeBX.resize(1);
  //  activeBX.coeffRef(0) = 0;
}



void run(std::string inputFile, std::string outFile)
{
  
  TFile *file2 = new TFile(inputFile.c_str());
  
 
  std::vector<double>* samples = new std::vector<double>;
  double amplitudeTruth;
  TTree *tree = (TTree*)file2->Get("Samples");
  tree->SetBranchAddress("amplitudeTruth",      &amplitudeTruth);
  tree->SetBranchAddress("samples",             &samples);
  int nentries = tree->GetEntries();
  
  
  
  
  TFile *fout;
  
  float time_shift = 13. ; //---- default is 13
  float pedestal_shift = 0.;
  
  
  float return_chi2 = -99;
  float best_pedestal = 0;
  float best_chi2 = 0;
  
  
  std::vector<TH1F*> v_pulses;
  std::vector<TH1F*> v_amplitudes_reco;
  
  std::cout << " outFile = " << outFile << std::endl;
  fout = new TFile(outFile.c_str(),"recreate");
  
  fout->cd();
  TTree* newtree = (TTree*) tree->CloneTree(0); //("RecoAndSim");
  newtree->SetName("RecoAndSim");
  
  std::vector <double> samplesReco;
  std::vector < std::vector<double> > complete_samplesReco;
  std::vector <double> complete_chi2;
  std::vector <double> complete_pedestal;
  
  int ipulseintime = 0;
  newtree->Branch("chi2",   &return_chi2, "chi2/F");
  newtree->Branch("samplesReco",   &samplesReco);
  newtree->Branch("ipulseintime",  ipulseintime,  "ipulseintime/I");
  
  newtree->Branch("complete_samplesReco",   &complete_samplesReco);
  newtree->Branch("complete_chi2",          &complete_chi2);
  newtree->Branch("complete_pedestal",          &complete_pedestal);
  newtree->Branch("best_pedestal",   &best_pedestal, "best_pedestal/F");
  newtree->Branch("best_chi2",   &best_chi2, "best_chi2/F");
  
  
  int totalNumberOfBxActive = 10;
  
  for (unsigned int ibx=0; ibx<totalNumberOfBxActive; ++ibx) {
    samplesReco.push_back(0.);
  }
  
  
  v_amplitudes_reco.clear();
  
  
  
  
  
  
  
  for(int ievt=0; ievt<nentries; ++ievt){
    tree->GetEntry(ievt);
    for(int i=0; i<NSAMPLES; i++){
      amplitudes[i] = samples->at(i);
    }
    
    double pedval = 0.;
    double pedrms = 1.0;
    PulseChiSqSNNLS pulsefunc;
    
    pulsefunc.disableErrorCalculation();
    bool status = pulsefunc.DoFit(amplitudes,noisecor,pedrms,activeBX,fullpulse,fullpulsecov);
    double chisq = pulsefunc.ChiSq();
    
    unsigned int ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<pulsefunc.BXs().rows(); ++ipulse) {
      if (pulsefunc.BXs().coeff(ipulse)==0) {
        ipulseintime = ipulse;
        break;
      }
    }
    double aMax = status ? pulsefunc.X()[ipulseintime] : 0.;
    //  double aErr = status ? pulsefunc.Errors()[ipulseintime] : 0.;
    
//     std::cout << " aMax = " << aMax << " amplitudeTruth = " << amplitudeTruth << "  chisq = " << chisq << std::endl;
    
    h01->Fill(aMax - amplitudeTruth);
  }
  std::cout << "  Mean of REC-MC = " << h01->GetMean() << " GeV" << std::endl;
  std::cout << "   RMS of REC-MC = " << h01->GetRMS() << " GeV" << std::endl;
}

void saveHist()
{
  
  fout->cd();
  h01->Write();
  fout->Close();
}



# ifndef __CINT__
int main(int argc, char** argv) {
  std::string inputFile = "data/samples_signal_10GeV_pu_0.root";
  if (argc>=2) {
    inputFile = argv[1];
  }
  
  std::string outFile = "output.root";
  if (argc>=3) {
    outFile = argv[2];
  }
  
  
  init();
  run(inputFile, outFile);
  saveHist();
  return 0;
}
# endif


