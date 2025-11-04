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

const int nTemplateBins = 9;
float pulseShapeTemplate[nTemplateBins];
std::vector<int> activeBXs = { -4, -3, -2, -1,  0,  1,  2 };

FullSampleVector fullpulse(FullSampleVector::Zero());
FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());
SampleMatrix noisecor(SampleMatrix::Zero());
BXVector activeBX;
SampleVector amplitudes(SampleVector::Zero());


void init()
{
  
  pSh.SetFNAMESHAPE("data/EmptyFileIdealPSphase2.root");
  pSh.Init();

  // intime sample is [3] // edm
  std::cout << "NFREQ = " << NFREQ << std::endl;
  for(int i=0; i<nTemplateBins; i++){
    
    //     double x = double( IDSTART + NFREQ * (i + 3) - WFLENGTH / 2);
    double x = double( NFREQ * i );    
    pulseShapeTemplate[i] = pSh.fShape(x);
    
  }
  //  for(int i=0; i<(NSAMPLES+2); i++) pulseShapeTemplate[i] /= pulseShapeTemplate[2];
  // 9 is the number of samples sufficient to cover the part non 0 of the pulse template
  // distance from min early BX (-4) to max late BX (+2) = 4*NFREQ + 16 + 2*NFREQ = 40 (NFREQ=4, fullpulse length)
  // shift from min early BX (-4) to first pulse sample (5) = 5 + 4 = 0
  for (int i=0; i<nTemplateBins; ++i) fullpulse(i+14) = pulseShapeTemplate[i];
    
  std::cout << " initialized fullpulse = " << std::endl << fullpulse << std::endl;
  
  for (int i=0; i<NSAMPLES; ++i) {
    for (int j=0; j<NSAMPLES; ++j) {
      int vidx = std::abs(j-i);
      noisecor(i,j) = pSh.corr(vidx);
    }
  }

  std::cout << " initialized noisecor = " << std::endl << noisecor << std::endl;

  activeBX.resize(activeBXs.size());
  for (unsigned int ibx=0; ibx<activeBX.size(); ++ibx) {
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
    
  std::cout << " outFile = " << outFile << std::endl;
  TFile *fout = new TFile(outFile.c_str(),"recreate");
  
  fout->cd();
  TTree* newtree = (TTree*) tree->CloneTree(0); //("RecoAndSim");
  newtree->SetName("RecoAndSim");
  
  float chisq;
  std::vector <double> samplesReco;
  
  int ipulseintime = 0;
  int nBins = nTemplateBins;
  newtree->Branch("chi2",   &chisq, "chi2/F");
  newtree->Branch("samplesReco",   &samplesReco);
  newtree->Branch("ipulseintime",  &ipulseintime,  "ipulseintime/I");
  newtree->Branch("activeBXs",     &activeBXs);
  newtree->Branch("nTemplateBins",   &nBins, "nTemplateBins/I");
  newtree->Branch("pulseShapeTemplate",   pulseShapeTemplate, "pulseShapeTemplate[nTemplateBins]/F");
  
  int totalNumberOfBxActive = activeBX.size();
  
  for (unsigned int ibx=0; ibx<totalNumberOfBxActive; ++ibx) {
    samplesReco.push_back(0.);
  }

  double pedval = 0.;
  double pedrms = 0.05;
  int maxshift = NPRESAMPLES + (*max_element(activeBX.begin(),activeBX.end())) * int(25./NFREQ);
  int minBX = *min_element(activeBX.begin(),activeBX.end());
  
  PulseChiSqSNNLS pulsefunc;
  pulsefunc.setNPresamples(NPRESAMPLES);
  pulsefunc.setNFREQ(NFREQ);
  pulsefunc.setMaxShift(maxshift);
  pulsefunc.disableErrorCalculation();

  
  for(int ievt=0; ievt<nentries; ++ievt){
    tree->GetEntry(ievt);
    for(int i=0; i<NSAMPLES; i++){
      amplitudes[i] = samples->at(i);
    }

    bool status = pulsefunc.DoFit(amplitudes,noisecor,pedrms,activeBX,fullpulse,fullpulsecov);
    chisq = pulsefunc.ChiSq();
    
    for (unsigned int ipulse=0; ipulse<pulsefunc.BXs().rows(); ++ipulse) {
      if (pulsefunc.BXs().coeff(ipulse)==0) {
        ipulseintime = ipulse;
        break;
      }
    }

    // std::cout << " ipulseintime = " << ipulseintime << std::endl;
    // std::cout << " pulsefunc.X() = " << std::endl << pulsefunc.X() << std::endl;
    // std::cout << " status = " << status << std::endl;
    // std::cout << " chi2 = " << chisq << std::endl;
    
    double aMax = status ? pulsefunc.X()[ipulseintime] : 0.;
    //  double aErr = status ? pulsefunc.Errors()[ipulseintime] : 0.;
    
    std::cout << " aMax = " << aMax << " amplitudeTruth = " << amplitudeTruth << "  chisq = " << chisq << std::endl;
    
    for (unsigned int ipulse=0; ipulse<pulsefunc.BXs().rows(); ++ipulse) {
      if (status) {
	samplesReco[ (int(pulsefunc.BXs().coeff(ipulse))) - minBX] = pulsefunc.X()[ ipulse ];
	// std::cout << "\t ipulse = " << ipulse << " idx = " << (int(pulsefunc.BXs().coeff(ipulse))) - minBX << "  ampli = " << pulsefunc.X()[ ipulse ] << std::endl;
      }
      else {
	samplesReco[ipulse] = -1;
      }
    }
    
    newtree->Fill();
    
  }

  newtree->Write();
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
  return 0;
}
# endif


