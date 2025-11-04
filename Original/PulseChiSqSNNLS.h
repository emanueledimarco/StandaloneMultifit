#ifndef PulseChiSqSNNLS_h
#define PulseChiSqSNNLS_h

#include "EigenMatrixTypes.h"
#include <set>
#include <array>

class PulseChiSqSNNLS {
public:
  
  typedef BXVector::Index Index;
  
  PulseChiSqSNNLS();
  ~PulseChiSqSNNLS();
  
  
  bool DoFit(const SampleVector &samples,
	     const SampleMatrix &samplecor,
	     double pederr,
	     const BXVector &bxs,
	     const FullSampleVector &fullpulse,
	     const FullSampleMatrix &fullpulsecov,
	     const SampleGainVector &gains = -1 * SampleGainVector::Ones(),
	     const SampleGainVector &badSamples = SampleGainVector::Zero());
  
  const SamplePulseMatrix &pulsemat() const { return _pulsemat; }
  const SampleMatrix &invcov() const { return _invcov; }
  
  const PulseVector &X() const { return _ampvecmin; }
  const PulseVector &Errors() const { return _errvec; }
  const BXVector &BXs() const { return _bxsmin; }
  
  double ChiSq() const { return _chisq; }
  void disableErrorCalculation() { _computeErrors = false; }
  void setNPresamples(int samples) { _npresamples = samples; }
  void setMaxShift(int maxshift) { _maxshift = maxshift; }
  void setNFREQ   ( float NFREQ )  { _NFREQ = NFREQ; }
  
protected:
  
  bool Minimize(const SampleMatrix &samplecor, double pederr, const FullSampleMatrix &fullpulsecov);
  bool NNLS();
  void NNLSUnconstrainParameter(Index idxp);
  void NNLSConstrainParameter(Index minratioidx);
  bool updateCov(const SampleMatrix &samplecor, double pederr, const FullSampleMatrix &fullpulsecov);
  double ComputeChiSq();
  double ComputeApproxUncertainty(unsigned int ipulse);
  
  
  SampleVector _sampvec;
  SampleMatrix _invcov;
  SamplePulseMatrix _pulsemat;
  PulseVector _ampvec;
  PulseVector _errvec;
  PulseVector _ampvecmin;
  
  SampleDecompLLT _covdecomp;
  
  BXVector _bxs;
  BXVector _bxsmin;
  unsigned int _npulsetot;
  unsigned int _nP;

  SamplePulseMatrix invcovp;
  PulseMatrix aTamat;
  PulseVector aTbvec;
  PulseVector updatework;

  PulseVector ampvecpermtest;
  
  double _chisq;
  bool _computeErrors;
  int _npresamples;
  int _maxshift;

  float _NFREQ;
  
};

#endif
