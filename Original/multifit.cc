//
// MultiFit amplitude reconstruction
// To run:
// > g++ -o Example06 Example06.cc PulseChiSqSNNLS.cc -std=c++11 `root-config --cflags --glibs`
// > ./Example06
//

#include <iostream>
#include <TString.h>
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
float templateCovariance[NSAMPLES][NSAMPLES];

std::vector<int> activeBXs = { -3, -2, -1,  0,  1,  2 };

FullSampleVector fullpulse(FullSampleVector::Zero());
FullSampleVector fullpulse_deriv(FullSampleVector::Zero());
FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());
SampleMatrix noisecor(SampleMatrix::Zero());
BXVector activeBX;
SampleVector amplitudes(SampleVector::Zero());
SampleGainVector gains(-1 * SampleGainVector::Ones());


void init()
{
  
  pSh.SetFNAMESHAPE("data/EmptyFileTestBeamPhase2.root");
  pSh.SetFNAMECOV("data/PulseCovarianceTestBeamPhase2.root");
  pSh.Init();

  // intime sample is [3] // edm
  for(int i=0; i<nTemplateBins; i++){
    
    //     double x = double( IDSTART + NFREQ * (i + 3) - WFLENGTH / 2);
    double x = double( NFREQ * i );
    pulseShapeTemplate[i] = pSh.fShape(x);
    
  }
  //  for(int i=0; i<(NSAMPLES+2); i++) pulseShapeTemplate[i] /= pulseShapeTemplate[2];
  // 9 is the number of samples sufficient to cover the part non 0 of the pulse template
  // distance from min early BX (-4) to max late BX (+2) = 4*NFREQ + 16 + 2*NFREQ = 40 (NFREQ=4, fullpulse length) (if min early BX =-3 then 3*4 + 16 + 2*4 = 36)
  // shift from min early BX (-4) to first pulse sample (5) = 5 + 4 = 0
  for (int i=0; i<nTemplateBins; ++i) fullpulse(i+14) = pulseShapeTemplate[i];

  for (int i = 0; i < nTemplateBins; ++i) {
    double x  = NFREQ * i;
    double dp = pSh.fShape(x + 0.5);
    double dm = pSh.fShape(x - 0.5);
    fullpulse_deriv(i + 14) = (dp - dm);
  }

  //std::cout << " initialized fullpulse = " << std::endl << fullpulse << std::endl;
  //std::cout << " initialized fullpulse_deriv = " << std::endl << fullpulse_deriv << std::endl;
  
  for(int i=0; i<NSAMPLES; i++) {
    for(int j=0; j<NSAMPLES; j++) {
      templateCovariance[i][j] = pSh.fCov(i,j);
      fullpulsecov(i + 14, j + 14) = templateCovariance[i][j];
    }
  }

  //  std::cout << " initialized fullpulsecov = " << std::endl << fullpulsecov << std::endl;

  for (int i=0; i<NSAMPLES; ++i) {
    for (int j=0; j<NSAMPLES; ++j) {
      int vidx = std::abs(j-i);
      noisecor(i,j) = pSh.corr(vidx);
    }
  }

  // std::cout << " initialized noisecor = " << std::endl << noisecor << std::endl;

  activeBX.resize(activeBXs.size());
  for (unsigned int ibx=0; ibx<activeBX.size(); ++ibx) {
    activeBX.coeffRef(ibx) = activeBXs[ibx];
  }
    //  activeBX.resize(1);
  //  activeBX.coeffRef(0) = 0;
}



void run(std::string inputFile, std::string outFile,
	 float sigmaNoiseScale=1,
	 int fitPedestal=0,
	 int maxEvents=-1)
{
  
  TFile *file2 = new TFile(inputFile.c_str());
 
  std::vector<double>* samples = new std::vector<double>;
  double amplitudeTruth;
  TTree *tree = (TTree*)file2->Get("Samples");
  tree->SetBranchAddress("amplitudeTruth",      &amplitudeTruth);
  tree->SetBranchAddress("samples",             &samples);
  int nentries = tree->GetEntries();
    
  TFile *fout = new TFile(outFile.c_str(),"recreate");
  
  fout->cd();
  TTree* newtree = (TTree*) tree->CloneTree(0); //("RecoAndSim");
  newtree->SetName("RecoAndSim");
  
  float chisq;
  std::vector <double> samplesReco;
  std::vector <double> pedestalsReco;
  
  int ipulseintime = 0;
  int nBins = nTemplateBins;
  newtree->Branch("chi2",   &chisq, "chi2/F");
  newtree->Branch("samplesReco",   &samplesReco);
  newtree->Branch("pedestalsReco",   &pedestalsReco);
  newtree->Branch("ipulseintime",  &ipulseintime,  "ipulseintime/I");
  newtree->Branch("activeBXs",     &activeBXs);
  newtree->Branch("nTemplateBins",   &nBins, "nTemplateBins/I");
  newtree->Branch("pulseShapeTemplate",   pulseShapeTemplate, "pulseShapeTemplate[nTemplateBins]/F");
  
  int totalNumberOfBxActive = activeBX.size();

  for (unsigned int ibx=0; ibx<totalNumberOfBxActive; ++ibx) {
    samplesReco.push_back(0.);
  }

  double pedval = 0.;
  double pedrms = sigmaNoiseScale*0.044;

  if (fitPedestal) gains = SampleGainVector::Zero(); // here decides n. pedestals to be fitted (1/gain)
  int ngains = gains.maxCoeff() + 1;
  for (int gainidx = 0; gainidx < ngains; ++gainidx) {
    SampleGainVector mask = gainidx * SampleGainVector::Ones();
    SampleVector pedestal = (gains.array() == mask.array()).cast<SampleVector::value_type>();
    if (pedestal.maxCoeff() > 0.) {
      pedestalsReco.push_back(0.);
    }
  }
  
  int maxshift = NPRESAMPLES + (*max_element(activeBX.begin(),activeBX.end())) * int(25./NFREQ);
  int minBX = *min_element(activeBX.begin(),activeBX.end());  
  
  // std::cout << "Total n of active BX (including pedestal) = " << totalNumberOfBxActive << std::endl;
  // std::cout << "Number of pedestals fitted = " << pedestalsReco.size() << std::endl;
  
  PulseChiSqSNNLS pulsefunc;
  pulsefunc.setNPresamples(NPRESAMPLES);
  pulsefunc.setNFREQ(NFREQ);
  pulsefunc.setMaxShift(maxshift);
  pulsefunc.disableErrorCalculation();

  
  for(int ievt=0; ievt<nentries; ++ievt){
    if (maxEvents>0 && ievt>=maxEvents) break;
    tree->GetEntry(ievt);
    for(int i=0; i<NSAMPLES; i++){
      amplitudes[i] = samples->at(i);
    }

    bool status = pulsefunc.DoFit(amplitudes,noisecor,pedrms,activeBX,fullpulse,fullpulse_deriv,fullpulsecov,gains);
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
    double aErr = status ? pulsefunc.Errors()[ipulseintime] : 0.;
    
    std::cout << " aMax = " << aMax << " amplitudeTruth = " << amplitudeTruth << "  chisq = " << chisq << std::endl;
    // std::cout << " aErr = " << aErr << std::endl;
    
    for (unsigned int ipulse=0; ipulse<pulsefunc.BXs().rows(); ++ipulse) {
      int iReco = (int(pulsefunc.BXs().coeff(ipulse)));
      if (status) {
        if (abs(iReco)<100) { 
          //          std::cout << "\t ipulse = " << ipulse << " idx = " << iReco << "  ampli = " << pulsefunc.X()[ ipulse ] << std::endl;
          samplesReco[iReco - minBX] = pulsefunc.X()[ ipulse ];
        } else if (iReco>=100) {
          std::cout << "pedestal[" << iReco-100 << "] = " << pulsefunc.X()[ ipulse ] << std::endl;
          pedestalsReco[iReco-100] = pulsefunc.X()[ ipulse ];
        } else {
          std::cout << " idx < 0 is for bad sample (e.g. slew rate). This should not happen, not turned on yet" << std::endl;
        }
      } else {
        if (iReco>=0 && iReco<100) {
          samplesReco[iReco - minBX] = -1;
        } else if (iReco>=100) {
          pedestalsReco[iReco-100] = -1;
        } else {
          std::cout << " idx < 0 is for bad sample (e.g. slew rate). This should not happen, not turned on yet" << std::endl;          
        }
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

  float sigmaNoiseScale = 1;
  if (argc>=3) {
    sigmaNoiseScale = atof(argv[2]);
  }

  int fitPedestal = 0;
  if (argc>=4) {
    fitPedestal = atoi(argv[3]);
  }

  char *suffix;
  if (argc>=5) {
    suffix = argv[4];
  } else {
    strcpy(suffix,"-");
  }

  int maxEvents = -1;
  if (argc>=6) {
    maxEvents = atoi(argv[5]);
  }
  
  std::cout << " input file = " << inputFile << std::endl;
  std::cout << " sigmaNoiseScale = "   << sigmaNoiseScale << std::endl;
  std::cout << " fitPedestal = " << fitPedestal << std::endl;
  std::cout << " suffix = " << suffix << std::endl;
  std::cout << " maxEvents = " << maxEvents << std::endl;

  
  TString outFile = Form("output_%s_noisescale%.2f_fitPed%d.root",suffix,sigmaNoiseScale,fitPedestal);

  std::cout << " output file = " << outFile.Data() << std::endl;
  
  init();
  run(inputFile, outFile.Data(),sigmaNoiseScale,fitPedestal,maxEvents);
  return 0;
}
# endif


