#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH2F.h>

#include <iostream>
#include <iomanip>

// #include "PulseParameters.h"



// total number of bunches in "LHC" bunch train
const int NBXTOTAL = 2800;

// length of a waveform in 1ns steps
// const int WFLENGTH  = 500;
// length of a waveform in 1ns/4 steps
const int WFLENGTH  = 500*4;

// number of samples per hit
const int NSAMPLES   = 16;

// number of pre-samples in impulse
int NPRESAMPLES = 6;

// distance between samples in 1ns steps
const float NFREQ      = 6.25;

// position of a 1st sample inside waveform
// const int IDSTART    = 180+13;
const int IDSTART    = 104;
// const int IDSTART    = 180;
// 

// CRRC shaping time in ns. It is used to calculate noise
// correlations. For QIE, set it to 1e-1
const double TAU = 43.0;




class Pulse{
  
  std::vector <double> _weights;
  std::vector <double> _mC;
  std::vector < std::vector< double > > _mL;
  
  //   double weights_[NSAMPLES];
  //   double mC_[NSAMPLES];
  //   double mL_[NSAMPLES][NSAMPLES];
  TGraph *_grPS;
  TH2F *_hCov;
  float _tMin;
  float _fPar0;
  float _fPar1;
  TFile *_filePS, *_fileCov;
  
  int _NSAMPLES;
  TString _FNAMESHAPE, _FNAMECOV;
  
  // distance between samples in 1ns steps
  //   a.k.a. ~1/sampling rate in ns   
  float _NFREQ; // = 25;
  
  // length of a waveform in 1ns steps
  int _WFLENGTH; // = 500;
  
  // position of a 1st sample inside waveform (in ns)
  float _IDSTART; //    = 180;
  
  // CRRC shaping time in ns. It is used to calculate noise
  // correlations. For QIE, set it to 1e-1
  float _TAU; // = 43.0;
  
  
public:
  
  Pulse();
  ~Pulse();
  
  void SetFNAMESHAPE ( std::string name );
  void SetFNAMECOV ( std::string name );
  void SetNSAMPLES ( int NSAMPLES );
  void SetNFREQ ( float NFREQ );
  void SetTAU ( float TAU );
  void SetWFLENGTH ( int WFLENGTH );
  void SetIDSTART ( float IDSTART );
  
  TGraph* grPS() {return _grPS; };
  float tMin() const { return _tMin; };
  float fPar0() const { return _fPar0; };
  float fPar1() const { return _fPar1; };
  float tau() const { return _TAU; };
  double weight(int i) const { return _weights.at(i); };
  double corr(int i) const { return _mC[i]; };
  double cholesky(int i, int j) const;
  
  void PrintNoise();
  void Init();
  void NoiseInit();
  void InitCholesky();
  void InitCorr();
  void SetNoiseCorrelationZero();
  void SetNoiseCorrelationMax();
  double fShape(double);
  float fCov(int i, int j);
  
};


Pulse::Pulse()
{
  //---- default
  SetFNAMESHAPE("data/EmptyFileCRRC43.root");
  SetFNAMECOV("data/PulseCovarianceTestBeamPhase2.root");
  SetNSAMPLES(16);
  SetNFREQ(6.25);
  SetIDSTART(104);
  SetTAU(43.0);
  SetWFLENGTH(208);
  
  _grPS = 0x0;
  _hCov = 0x0;
}


Pulse::~Pulse()
{
}



double Pulse::cholesky ( int i, int j ) const {
  //  std::cout << " i,j = " << i << " , " << j << std::endl;
  //  std::cout << " _mL.size() = " << _mL.size() << std::endl;
  //  std::cout << " _mL.at(i).size() = " << _mL.at(i).size() << std::endl;
  //  std::cout << " mL = " << _mL.at(i).at(j) << std::endl;
  return _mL.at(i).at(j);
}


void Pulse::SetIDSTART ( float IDSTART ) {
  _IDSTART = IDSTART;
}


void Pulse::SetTAU ( float TAU ) {
  _TAU = TAU;
}


void Pulse::SetWFLENGTH ( int WFLENGTH ) {
  _WFLENGTH = WFLENGTH;
}


void Pulse::SetNFREQ ( float NFREQ ) {
  _NFREQ = NFREQ;
}


void Pulse::SetNSAMPLES ( int NSAMPLES ) {
  _NSAMPLES = NSAMPLES;
}

void Pulse::SetFNAMESHAPE ( std::string name ) {
  _FNAMESHAPE = Form ("%s", name.c_str());
}

void Pulse::SetFNAMECOV ( std::string name ) {
  _FNAMECOV = Form ("%s", name.c_str());
}


void Pulse::Init() {
  
  _filePS = new TFile(_FNAMESHAPE.Data());
  _grPS = (TGraph*) ((TGraph*)_filePS->Get("PulseShape/grPulseShape")) -> Clone();

  _fileCov = new TFile(_FNAMECOV.Data());
  _hCov = (TH2F*) ((TH2F*)_fileCov->Get("PulseCovariance")) -> Clone();
  
  // In-time sample is i=5
  for(int i=0; i<_NSAMPLES; i++){
    double x = double( _IDSTART + _NFREQ * i - _WFLENGTH / 2);
    _weights.push_back( fShape(x) );
  }
  
  NoiseInit(); 
}



double Pulse::fShape(double x) {
  
  if ( _grPS !=0 && x > 0.) {
      return _grPS->Eval(x);
  }
  else {
    return 0.;
  }
}

float Pulse::fCov(int i, int j) {
  return _hCov->GetBinContent(i+1,j+1);
  //if ( i == j ) return 0.;
  //if ( i>=0 && i<_hCov->GetNbinsX()) return _hCov->GetBinContent(i+1,j+1);
  //return 0.;
}

void Pulse::SetNoiseCorrelationZero() {
  _mC.clear();
  _mC.push_back(1.0);
  for(int i=1; i<_NSAMPLES; i++){
    _mC.push_back(0.0);
  }
  
  _mL.clear();
  // initialize
  for(int i=0; i<_NSAMPLES; ++i){
    std::vector<double> temp_mL;
    for(int j=0; j<_NSAMPLES; ++j){
      temp_mL.push_back(0);
    }
    _mL.push_back(temp_mL);
  }
  for(int i=0; i<_NSAMPLES; i++){
    _mL.at(i).at(i) = 1;
  }
}

void Pulse::SetNoiseCorrelationMax() {
  _mC.clear();
  for(int i=0; i<_NSAMPLES; i++){
    _mC.push_back(1.0);
  }
  
  _mL.clear();
  // initialize
  for(int i=0; i<_NSAMPLES; ++i){
    std::vector<double> temp_mL;
    for(int j=0; j<_NSAMPLES; ++j){
      temp_mL.push_back(0);
    }
    _mL.push_back(temp_mL);
  }
  for(int i=0; i<_NSAMPLES; i++){
    _mL.at(i).at(0) = 1;
  }
}


void Pulse::InitCorr() {
  
  //  std::cout << " >> Pulse::NoiseInit " << std::endl;
  _mC.clear();
  
  for(int i=0; i<_NSAMPLES; i++){
    double y = 1. - exp( -double(_NFREQ * i) / (sqrt(2.) * _TAU));
    _mC.push_back( 1. - y * y);
//     _mC.push_back( 1. );
//     _mC.push_back( 0.999 );
//        _mC.push_back( 0.0 );
  }
}


void Pulse::InitCholesky() {
  
  _mL.clear();
  
  // initialize
  for(int i=0; i<_NSAMPLES; ++i){
    std::vector<double> temp_mL;
    for(int j=0; j<_NSAMPLES; ++j){
      temp_mL.push_back(0);
    }
    _mL.push_back(temp_mL);
  }
  
  // decomposition
  if (_mC.at(0)< 0) {
    //    std::cout << "  ----> _mC.at(0) = " << _mC.at(0) << std::endl;
  }
  _mL.at(0).at(0) = sqrt(_mC.at(0));
  for( int col=1; col<_NSAMPLES; col++){
    _mL.at(0).at(col) = 0;
  }
  for( int row=1; row<_NSAMPLES; row++){
    for( int col=0; col<row; col++ ){
      double sum1 = 0;
      int m=abs(row-col);
      for( int k=0; k<col; ++k) sum1 += _mL.at(row).at(k)*_mL.at(col).at(k);
      _mL.at(row).at(col) = (_mC.at(m) - sum1)/_mL.at(col).at(col);
    }
    double sum2 = 0;
    for( int k=0; k<row; ++k) sum2 += _mL.at(row).at(k)*_mL.at(row).at(k);
    //     std::cout << "  ----> _mC.at(0) - sum2 = " << _mC.at(0) << " - " << sum2 << " = " << _mC.at(0) - sum2 << std::endl;
    
    _mL.at(row).at(row) = sqrt( _mC.at(0) - sum2 );
    for( int col=row+1; col<_NSAMPLES; col++ ) _mL.at(row).at(col) = 0;
  }
}





void Pulse::PrintNoise () {
  
  for (int i=0; i < _mL.size(); i++ ) {
    for (int j=0; j < _mL.at(i).size(); j++ ) {
      std::cout << std::setw(5) << std::setprecision(3) << _mL.at(i).at(j) << "   " ;
    }
    std::cout << std::endl;
  } 
}  


void Pulse::NoiseInit() {
  InitCorr();
  InitCholesky();
}


