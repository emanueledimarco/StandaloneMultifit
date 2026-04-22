//
// MultiFit amplitude reconstruction
// To run:
// > g++ -o Example06 Example06.cc PulseChiSqSNNLS.cc -std=c++11 `root-config --cflags --glibs`
// > ./Example06
//

#include <iostream>
#include <TString.h>
#include "PulseChiSqSNNLS.h"

#include "TTree.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TFile.h"

using namespace std;

 Pulse pSh;

const int nTemplateBins = 11;
float templateCovariance[NSAMPLES][NSAMPLES];
float pulseShapeTemplate[nTemplateBins], signalTemplateError[nTemplateBins];


std::vector<int> activeBXs = { -3, -2, -1,  0,  1, 2};

FullSampleVector fullpulse(FullSampleVector::Zero());
FullSampleVector fullpulse_deriv(FullSampleVector::Zero());
FullSampleVector fullpulse_deriv2(FullSampleVector::Zero());
FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());
SampleMatrix noisecor(SampleMatrix::Zero());
BXVector activeBX;
SampleVector amplitudes(SampleVector::Zero());
SampleGainVector gains(-1 * SampleGainVector::Ones());

FullSampleVector fullpulse_signal_template_error(FullSampleVector::Zero());


void init()
{

  pSh.SetFNAMESHAPE("data/TestBeamPhase2_PS_coeffs.txt");
  pSh.SetFNAMECOV("data/PulseCovarianceTestBeamPhase_withFlatSignalError.root");
  pSh.Init();

  pSh.SetNoiseCorrelationZero();

  std::cout << std::endl << std::endl << "printing read PS" << std::endl;
  //   intime sample is [3] // edm
  for(int i=0; i<nTemplateBins; i++){
    //     double x = double( IDSTART + NFREQ * (i + 3) - WFLENGTH / 2);
    double x = double( NFREQ * i - PULSESHAPE_SHIFT );
    pulseShapeTemplate[i] = pSh.fShape(x);
    std::cout << pulseShapeTemplate[i] << std::endl;
    signalTemplateError[i] = pSh.fSignalShapeError(x);
  }

  for (int i=0; i<nTemplateBins; ++i){
    fullpulse(i+14) = pulseShapeTemplate[i];
    fullpulse_signal_template_error(i+14) = signalTemplateError[i];
  }

  for (int i = 0; i < nTemplateBins; ++i) {
    double x  = NFREQ * i;
    double dp = pSh.fShape(x + 0.1 - PULSESHAPE_SHIFT );
    double dm = pSh.fShape(x - 0.1 - PULSESHAPE_SHIFT );
    fullpulse_deriv(i + 14) = (dp - dm)/0.2;
    fullpulse_deriv2(i + 14) = (dp + dm - 2*pulseShapeTemplate[i])/0.01;
  }

  std::cout << " initialized fullpulse = " << std::endl << fullpulse << std::endl;
  std::cout << " initialized fullpulse_deriv = " << std::endl << fullpulse_deriv << std::endl;

  for(int i=0; i<NSAMPLES; i++) {
    for(int j=0; j<NSAMPLES; j++) {
      templateCovariance[i][j] = pSh.fCov(i,j);
      fullpulsecov(i + 14, j + 14) = templateCovariance[i][j];
    }
  }


  std::cout << " initialized fullpulsecov = " << std::endl << fullpulsecov << std::endl;

  for (int i=0; i<NSAMPLES; ++i) {
    for (int j=0; j<NSAMPLES; ++j) {
      int vidx = std::abs(j-i);
      noisecor(i,j) = pSh.corr(vidx);
      //if ( (i == 9) && (j==9) ) noisecor(i,j) = 100000.;
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
  std::vector <double> samplesReco, timeReco;
  std::vector <double> pedestalsReco, normResVector, absResVector;

  int ipulseintime = 0;
  int nBins = nTemplateBins;
  newtree->Branch("chi2",   &chisq, "chi2/F");
  newtree->Branch("samplesReco",   &samplesReco);
  newtree->Branch("timeReco",   &timeReco);
  newtree->Branch("pedestalsReco",   &pedestalsReco);
  newtree->Branch("absResVector",   &absResVector);
  newtree->Branch("normResVector",   &normResVector);
  newtree->Branch("ipulseintime",  &ipulseintime,  "ipulseintime/I");
  newtree->Branch("activeBXs",     &activeBXs);
  newtree->Branch("nTemplateBins",   &nBins, "nTemplateBins/I");
  newtree->Branch("pulseShapeTemplate",   pulseShapeTemplate, "pulseShapeTemplate[nTemplateBins]/F");
  
  int totalNumberOfBxActive = activeBX.size();

  for (unsigned int ibx=0; ibx<totalNumberOfBxActive; ++ibx) {
    samplesReco.push_back(0.);
    timeReco.push_back(0.);
  }

  for(int i=0; i<NSAMPLES; i++){
     absResVector.push_back(0.);
     normResVector.push_back(0.);
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

  SampleGainVector badsamples = SampleGainVector::Zero();
  //badsamples[9] = 1;

  for(int ievt=0; ievt<nentries; ++ievt){
    if (maxEvents>0 && ievt>=maxEvents) break;
    tree->GetEntry(ievt);
    for(int i=0; i<NSAMPLES; i++){
      amplitudes[i] = samples->at(i) + (fitPedestal ? FIXED_PEDESTAL : 0);
    }

    bool status = pulsefunc.DoFit(amplitudes,noisecor,pedrms,activeBX,fullpulse,fullpulse_deriv,fullpulsecov,fullpulse_signal_template_error,pSh,gains,badsamples);
    chisq = pulsefunc.ChiSq();

    SampleVector normResVec = pulsefunc.NormRes();
    SampleVector absResVec = pulsefunc.AbsRes();

    for(int i=0; i<NSAMPLES; i++){
      normResVector[i] = normResVec(i, 0);
      absResVector[i] = absResVec(i, 0);
   }

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
    double time = status ? pulsefunc.T()[ipulseintime] : 0.;
    
    //std::cout << " aMax = " << aMax << " amplitudeTruth = " << amplitudeTruth << "  chisq = " << chisq << std::endl;
    // std::cout << " aErr = " << aErr << std::endl;
    double pedestal_offset = FIXED_PEDESTAL;
    for (unsigned int ipulse=0; ipulse<pulsefunc.BXs().rows(); ++ipulse) {
      int iReco = (int(pulsefunc.BXs().coeff(ipulse)));
      if (status) {
        if (abs(iReco)<PEDESTAL_BX_OFFSET) {
          //          std::cout << "\t ipulse = " << ipulse << " idx = " << iReco << "  ampli = " << pulsefunc.X()[ ipulse ] << std::endl;
          samplesReco[iReco - minBX] = pulsefunc.X()[ ipulse ];
          timeReco[iReco - minBX] = pulsefunc.T()[ ipulse ];
        } else if (iReco>=PEDESTAL_BX_OFFSET) {
          pedestalsReco[iReco-PEDESTAL_BX_OFFSET] = pulsefunc.X()[ ipulse ] - pedestal_offset;
          pedestal_offset = 0;
          std::cout << "pedestal[" << iReco-PEDESTAL_BX_OFFSET << "] = " << pulsefunc.X()[ ipulse ] << std::endl;
        } else {
          std::cout << " idx < 0 is for bad sample (e.g. slew rate). This should not happen, not turned on yet" << std::endl;
        }
      } else {
        if (iReco>=0 && iReco<PEDESTAL_BX_OFFSET) {
          samplesReco[iReco - minBX] = -1;
          timeReco[iReco - minBX] = -1;
        } else if (iReco>=PEDESTAL_BX_OFFSET) {
          pedestalsReco[iReco-PEDESTAL_BX_OFFSET] = -1;
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


