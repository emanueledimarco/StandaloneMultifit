//#define PulseChiSqSNNLS_cxx
#include "PulseChiSqSNNLS.h"
#include <math.h>
#include <iostream>

PulseChiSqSNNLS::PulseChiSqSNNLS() :
_chisq(0.),
_deltachisq(0.),
_computeErrors(true),
_npresamples(3),
_maxshift(7)
{
  
  Eigen::initParallel();
  _invcov.setZero();
}  

PulseChiSqSNNLS::~PulseChiSqSNNLS() {
  
}

void eigen_solve_submatrix(PulseMatrix &mat, PulseVector &invec, PulseVector &outvec, unsigned NP) {
  using namespace Eigen;
  switch (NP) {  // pulse matrix is always square.
    case 10: {
      Matrix<double, 10, 10> temp = mat.topLeftCorner<10, 10>();
      outvec.head<10>() = temp.ldlt().solve(invec.head<10>());
    } break;
    case 9: {
      Matrix<double, 9, 9> temp = mat.topLeftCorner<9, 9>();
      outvec.head<9>() = temp.ldlt().solve(invec.head<9>());
    } break;
    case 8: {
      Matrix<double, 8, 8> temp = mat.topLeftCorner<8, 8>();
      outvec.head<8>() = temp.ldlt().solve(invec.head<8>());
    } break;
    case 7: {
      Matrix<double, 7, 7> temp = mat.topLeftCorner<7, 7>();
      outvec.head<7>() = temp.ldlt().solve(invec.head<7>());
    } break;
    case 6: {
      Matrix<double, 6, 6> temp = mat.topLeftCorner<6, 6>();
      outvec.head<6>() = temp.ldlt().solve(invec.head<6>());
    } break;
    case 5: {
      Matrix<double, 5, 5> temp = mat.topLeftCorner<5, 5>();
      outvec.head<5>() = temp.ldlt().solve(invec.head<5>());
    } break;
    case 4: {
      Matrix<double, 4, 4> temp = mat.topLeftCorner<4, 4>();
      outvec.head<4>() = temp.ldlt().solve(invec.head<4>());
    } break;
    case 3: {
      Matrix<double, 3, 3> temp = mat.topLeftCorner<3, 3>();
      outvec.head<3>() = temp.ldlt().solve(invec.head<3>());
    } break;
    case 2: {
      Matrix<double, 2, 2> temp = mat.topLeftCorner<2, 2>();
      outvec.head<2>() = temp.ldlt().solve(invec.head<2>());
    } break;
    case 1: {
      Matrix<double, 1, 1> temp = mat.topLeftCorner<1, 1>();
      outvec.head<1>() = temp.ldlt().solve(invec.head<1>());
    } break;
    default:
      throw std::domain_error("MultFitWeirdState: Weird number of pulses encountered in multifit, module is configured incorrectly!");
  }
}

bool PulseChiSqSNNLS::DoFit(const SampleVector &samples,
			    const SampleMatrix &samplecor,
			    double pederr,
			    const BXVector &bxs,
			    const FullSampleVector &fullpulse,
          const FullSampleVector &fullpulse_deriv,
			    const FullSampleMatrix &fullpulsecov,
			    const Pulse &pSh,
          const SampleGainVector &gains,
          const SampleGainVector &badSamples
			    ) {

  _pSh = pSh;

  const unsigned int nsample = SampleVector::RowsAtCompileTime;
  const unsigned int npulse = bxs.rows();

  _sampvec = samples;

  _normResVec = SampleVector::Zero();
  _absResVec = SampleVector::Zero();

  _bxs = bxs;
  _pulsemat.resize(Eigen::NoChange, npulse);
  _pulsemat_t.resize(Eigen::NoChange, npulse);

  //construct dynamic pedestals if applicable
  int ngains = gains.maxCoeff() + 1;
  int nPedestals = 0;
  for (int gainidx = 0; gainidx < ngains; ++gainidx) {
    SampleGainVector mask = gainidx * SampleGainVector::Ones();
    SampleVector pedestal = (gains.array() == mask.array()).cast<SampleVector::value_type>();
    if (pedestal.maxCoeff() > 0.) {
      ++nPedestals;
      _bxs.resize(npulse + nPedestals);
      _bxs[npulse + nPedestals - 1] = 100 + gainidx;  //bx values >=100 indicate dynamic pedestals
      _pulsemat.resize(Eigen::NoChange, npulse + nPedestals);
       _pulsemat_t.resize(Eigen::NoChange, npulse + nPedestals);
      _pulsemat.col(npulse + nPedestals - 1) = pedestal;
    }
  }

  //construct negative step functions for saturated or potentially slew-rate-limited samples
  for (int isample = 0; isample < SampleVector::RowsAtCompileTime; ++isample) {
    if (badSamples.coeff(isample) > 0) {
      SampleVector step = SampleVector::Zero();
      //step correction has negative sign for saturated or slew-limited samples which have been forced to zero
      step[isample] = -1.;

      ++nPedestals;
      _bxs.resize(npulse + nPedestals);
      _bxs[npulse + nPedestals - 1] =
          -100 - isample;  //bx values <=-100 indicate step corrections for saturated or slew-limited samples
      _pulsemat.resize(Eigen::NoChange, npulse + nPedestals);
      _pulsemat_t.resize(Eigen::NoChange, npulse + nPedestals);
      _pulsemat.col(npulse + nPedestals - 1) = step;
    }
  }

  _npulsetot = npulse + nPedestals;
  
  _ampvec = PulseVector::Zero(_npulsetot);
  _errvec = PulseVector::Zero(_npulsetot);
  _nP = 0;
  _chisq = 0.;
  _deltachisq = 0.;
  
  aTamat.resize(_npulsetot, _npulsetot);

  _time      = PulseVector::Zero(_npulsetot);  // Δt = 0 for all pulses initially
  _timeErr   = PulseVector::Zero(_npulsetot);
  _timeActive = Eigen::VectorXi::Zero(_npulsetot); // default all inactive

  //initialize pulse template matrix
  for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
    int bx = _bxs.coeff(ipulse);
    int firstsamplet = std::max(0,bx * int(25./_NFREQ) + _npresamples);
    int offset = _maxshift - _npresamples - bx*int(25./_NFREQ);
    _pulsemat.col(ipulse) = fullpulse.segment<SampleVector::RowsAtCompileTime>(offset);
    _pulsemat_t.col(ipulse) = fullpulse_deriv.segment<SampleVector::RowsAtCompileTime>(offset);

  }

  for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
    //std::cout << "pulsemat at ipulse = " << ipulse << ": " << _pulsemat.col(ipulse) << std::endl;
    //std::cout << "pulsemat_t at ipulse = " << ipulse << ": " << _pulsemat_t.col(ipulse) << std::endl;
  }

  // std::cout << "pulsemat after pedestal + bad samples additions = " << std::endl << _pulsemat << std::endl;

  //unconstrain pedestals already for first iteration since they should always be non-zero
  if (nPedestals > 0) {
    for (int i = 0; i < _bxs.rows(); ++i) {
      int bx = _bxs.coeff(i);
      if (bx >= 100) {
        NNLSUnconstrainParameter(i);
      }
    }
  }

  //do the actual fit
  bool status = Minimize(samplecor,pederr,fullpulsecov);
  _ampvecmin = _ampvec;
  
  // std::cout << " _sampvec = " << _sampvec << std::endl;
  // std::cout << " bxs = " << bxs << std::endl;
  // std::cout << " fullpulse = " << fullpulse << std::endl;
  // std::cout << " _ampvecmin = " << _ampvecmin << std::endl;
  
  _bxsmin = _bxs;
  
  if (!status) return status;
  
  //   std::cout << " _computeErrors = " << _computeErrors << std::endl;
  
  if(!_computeErrors) return status;
  
  //compute MINOS-like uncertainties for in-time amplitude
  bool foundintime = false;
  unsigned int ipulseintime = 0;
  //   std::cout << " npulse = " << npulse << std::endl;
  for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
    //     std::cout << " _bxs.coeff( " << ipulse << "::" << npulse << " ) = " << _bxs.coeff(ipulse) << std::endl;
    if (_bxs.coeff(ipulse)==0) {
      ipulseintime = ipulse;
      foundintime = true;
      break;
    }
  }
  // std::cout << " foundintime = " << foundintime << std::endl;
  if (!foundintime) return status;
  
  
  
  const unsigned int ipulseintimemin = ipulseintime;
  
  double approxerr = ComputeApproxUncertainty(ipulseintime);
  double chisq0 = _chisq;
  double x0 = _ampvecmin[ipulseintime];
  
  //move in time pulse first to active set if necessary
  if (ipulseintime<_nP) {
    _pulsemat.col(_nP-1).swap(_pulsemat.col(ipulseintime));
    _pulsemat_t.col(_nP-1).swap(_pulsemat_t.col(ipulseintime));
    std::swap(_ampvec.coeffRef(_nP-1),_ampvec.coeffRef(ipulseintime));
    std::swap(_bxs.coeffRef(_nP-1),_bxs.coeffRef(ipulseintime));
    ipulseintime = _nP - 1;
    --_nP;    
  }
  
  
  SampleVector pulseintime = _pulsemat.col(ipulseintime);
  _pulsemat.col(ipulseintime).setZero();
  
  //two point interpolation for upper uncertainty when amplitude is away from boundary
  double xplus100 = x0 + approxerr;
  _ampvec.coeffRef(ipulseintime) = xplus100;
  _sampvec = samples - _ampvec.coeff(ipulseintime)*pulseintime;

  status &= Minimize(samplecor,pederr,fullpulsecov);
  if (!status) return status;
  double chisqplus100 = ComputeChiSq();
  
  double sigmaplus = std::abs(xplus100-x0)/sqrt(chisqplus100-chisq0);
  
  //if amplitude is sufficiently far from the boundary, compute also the lower uncertainty and average them
  if ( (x0/sigmaplus) > 0.5 ) {
    for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
      if (_bxs.coeff(ipulse)==0) {
        ipulseintime = ipulse;
        break;
      }
    }    
    double xminus100 = std::max(0.,x0-approxerr);
    _ampvec.coeffRef(ipulseintime) = xminus100;
    _sampvec = samples - _ampvec.coeff(ipulseintime)*pulseintime;
    status &= Minimize(samplecor,pederr,fullpulsecov);
    if (!status) return status;
    double chisqminus100 = ComputeChiSq();
    
    double sigmaminus = std::abs(xminus100-x0)/sqrt(chisqminus100-chisq0);
    _errvec[ipulseintimemin] = 0.5*(sigmaplus + sigmaminus);
    
  }
  else {
    _errvec[ipulseintimemin] = sigmaplus;
  }
  
  _chisq = chisq0;  
  
  return status;
  
}

bool PulseChiSqSNNLS::Minimize(const SampleMatrix &samplecor, double pederr, const FullSampleMatrix &fullpulsecov) {
  const unsigned int npulse = _bxs.rows();
  
  const int maxiter = 50;
  int iter = 0;
  bool status = false;
  while (true) {    
    
    if (iter>=maxiter) {
      std::cout << "PulseChiSqSNNLS::Minimize ===> " << "Max Iterations reached at iter " << iter <<  std::endl;
      std::cout << " maxiter =  " << iter << " :: " << maxiter << std::endl;
      break;
    }    
    
    status = updateCov(samplecor,pederr,fullpulsecov);    
    if (!status) break;    
    if (npulse > 1) {
      status = NNLS();
    } else {
      //special case for one pulse fit (performance optimized)
      status = OnePulseMinimize();
    }
    if (!status) break;
    
    double chisqnow = ComputeChiSq();
    double deltachisq = chisqnow-_chisq;

    // std::cout << "Iter = " << iter << "  chisq now = " << chisqnow <<  "   deltachisq = " << std::abs(deltachisq) << std::endl;
    // std::cout << "N active pulses = " << _nP << std::endl;
    
    _chisq = chisqnow;
    if (std::abs(deltachisq)<1e-3) {
      break;
    }

    // std::cout << "              _deltachisq = " << _deltachisq << " doubleDelta = " << std::abs(std::abs(_deltachisq)-std::abs(deltachisq)) << std::endl;
    // to avoid bouncing between two degenerate solutions
    if (std::abs(std::abs(_deltachisq)-std::abs(deltachisq))<1e-3) {
      break;
    }
    _deltachisq = deltachisq;
    
    ++iter;    
  }  


  TimingSignalRefit();
//
//  const int nTemplateBins = 9;
//  float pulseShapeTemplate[nTemplateBins];
//  FullSampleVector fullpulse(FullSampleVector::Zero());
//  FullSampleVector fullpulse_deriv(FullSampleVector::Zero());
//
//  // intime sample is [3] // edm
//  for(int i=0; i<nTemplateBins; i++){
//    //     double x = double( IDSTART + NFREQ * (i + 3) - WFLENGTH / 2);
//    double x = double( NFREQ * i - PULSESHAPE_SHIFT );
//    pulseShapeTemplate[i] = _pSh.fShape(x);
//  }
//
//  for (int i=0; i<nTemplateBins; ++i) fullpulse(i+14) = pulseShapeTemplate[i];
//
//  for (int i = 0; i < nTemplateBins; ++i) {
//    double x  = NFREQ * i;
//    double dp = _pSh.fShape(x + 0.1 - PULSESHAPE_SHIFT );
//    double dm = _pSh.fShape(x - 0.1 - PULSESHAPE_SHIFT );
//    fullpulse_deriv(i + 14) = (dp - dm)/0.2;
//  }
//
//
//  //initialize pulse template matrix
//  for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
//    int bx = _bxs.coeff(ipulse);
//    int firstsamplet = std::max(0,bx * int(25./_NFREQ) + _npresamples);
//    int offset = _maxshift - _npresamples - bx*int(25./_NFREQ);
//    _pulsemat.col(ipulse) = fullpulse.segment<SampleVector::RowsAtCompileTime>(offset);
//    _pulsemat_t.col(ipulse) = fullpulse_deriv.segment<SampleVector::RowsAtCompileTime>(offset);
//
//  }
//

  //initialize pulse template matrix
  for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
     std::cout << "_time[ipulse]: " << _time[ipulse] << std::endl;
    _pulsemat.col(ipulse) -= _time[ipulse]*_pulsemat_t.col(ipulse);
  }

  TimingSignalRefit();

  return status;

}

void PulseChiSqSNNLS::TimingSignalRefit() {

    // -------------------------
    // Identify in-time pulse
    // -------------------------
    unsigned int ipulseSignal = 0;
    bool foundSignal = false;
    for (unsigned int ip = 0; ip < _bxs.rows(); ++ip) {
        if (_bxs.coeff(ip) == 0) {
            ipulseSignal = ip;
            foundSignal = true;
            break;
        }
    }
    if (!foundSignal) return;

    // -------------------------
    // Build residual including all other pulses
    // -------------------------
    SampleVector res = _sampvec;
    for (unsigned int ip = 0; ip < _bxs.rows(); ++ip) {
        if (ip != ipulseSignal) res -= _ampvec.coeff(ip) * _pulsemat.col(ip);
    }

    // Extract pulse and derivative for the in-time pulse
    Eigen::VectorXd P  = _pulsemat.col(ipulseSignal);
    Eigen::VectorXd Pd = _pulsemat_t.col(ipulseSignal);

    // Cholesky of the covariance
    Eigen::MatrixXd L = _covdecomp.matrixL();

    // Whiten vectors
    Eigen::VectorXd Pw  = L.triangularView<Eigen::Lower>().solve(P);
    Eigen::VectorXd Pdw = L.triangularView<Eigen::Lower>().solve(Pd);

    // Compute weighted cross-term
    double cross = Pw.dot(Pdw);
    double corr  = cross / (Pw.norm() * Pdw.norm());

    std::cout << "[Diagnostics] Weighted cross-term = " << cross
              << ", correlation = " << corr << std::endl;


    // -------------------------
    // Design matrix for active fit: columns = [P, P']
    // -------------------------
    Eigen::Matrix<double, SampleVector::RowsAtCompileTime, 2> M;
    M.col(0) = _pulsemat.col(ipulseSignal);
    M.col(1) = _pulsemat_t.col(ipulseSignal);

    // -------------------------
    // NNLS 2D solve
    // -------------------------
    Eigen::Vector2d x = Eigen::Vector2d::Zero();       // solution [A, C]
    bool activeA = false;
    bool activeC_upper = false;
    bool activeC_lower = false;

    const int maxIter = 10;
    for (int iter = 0; iter < maxIter; ++iter) {

        // Solve unconstrained LS on free variables
        Eigen::MatrixXd M_free( SampleVector::RowsAtCompileTime, 2 );
        Eigen::Vector2d x_free;

        int nFree = 0;
        Eigen::MatrixXd M_sub( SampleVector::RowsAtCompileTime, 0 );
        Eigen::Vector2d rhs_sub;
        std::vector<int> freeIdx;

        if (!activeA) freeIdx.push_back(0);
        if (!activeC_upper && !activeC_lower) freeIdx.push_back(1);

        nFree = int(freeIdx.size());
        if (nFree == 0) break; // all active, solution fixed

        M_sub.resize( SampleVector::RowsAtCompileTime, nFree );
        for (int j=0; j<nFree; ++j) M_sub.col(j) = M.col(freeIdx[j]);

        // Solve unconstrained on free subset
        Eigen::VectorXd rhs = _covdecomp.matrixL().solve(res);
        Eigen::VectorXd M_rhs = M_sub.transpose() * rhs;
        Eigen::MatrixXd M_mat = M_sub.transpose() * _covdecomp.matrixL().solve(M_sub);

        Eigen::VectorXd xsol = M_mat.ldlt().solve(M_rhs);

        // Assign back to x
        for (int j=0; j<nFree; ++j) x(freeIdx[j]) = xsol(j);

        // -------------------------
        // Check constraints
        // -------------------------
        bool violated = false;

        if (x(0) < 0.) { x(0) = 0.; activeA = true; violated = true; }
        double Cmin = -3.0 * x(0);
        double Cmax =  3.0 * x(0);
        if (x(1) < Cmin) { x(1) = Cmin; activeC_lower = true; violated = true; }
        if (x(1) > Cmax) { x(1) = Cmax; activeC_upper = true; violated = true; }

        if (!violated) break; // solution feasible
    }

    // -------------------------
    // Store results
    // -------------------------
    _ampvec.coeffRef(ipulseSignal) = x(0);
    _time.coeffRef(ipulseSignal) = - x(1) / x(0); // Δt = -C / A
    _timeActive(ipulseSignal) = 1;

    // update χ²
    _chisq = ComputeChiSq();

    std::cout << "TimingSignalRefit -> A = " << _ampvec[ipulseSignal]
              << ", Δt = " << _time[ipulseSignal] << std::endl;
}



bool PulseChiSqSNNLS::updateCov(const SampleMatrix &samplecor, double pederr, const FullSampleMatrix &fullpulsecov) {
  
  // std::cout << " updateCov " << std::endl;
  
  const unsigned int nsample = SampleVector::RowsAtCompileTime;
  const unsigned int npulse = _bxs.rows();

  // std::cout << "_bxs = " << std::endl << _bxs << "   npulse = " << npulse << std::endl;

  _invcov.triangularView<Eigen::Lower>() = (pederr*pederr)*samplecor;

  //std::cout << " invcov (only noise) = " << std::endl << _invcov << std::endl;

  for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
    if (_ampvec.coeff(ipulse)==0.) continue;
    int bx = _bxs.coeff(ipulse);
    if (std::abs(bx) >= 100)
      continue;  //no contribution to covariance from pedestal or saturation/slew step correction

    int firstsamplet = std::max(0,bx * int(25./_NFREQ) + _npresamples);
    int offset = _maxshift - _npresamples - bx*int(25./_NFREQ);

    double ampsq = _ampvec.coeff(ipulse)*_ampvec.coeff(ipulse);
    // std::cout << "     >>> ipulse = " << ipulse << "    ampsq = " << ampsq << std::endl;
    const unsigned int nsamplepulse = nsample-firstsamplet;

    const auto &block_allocated = fullpulsecov.block(firstsamplet+offset,firstsamplet+offset,nsamplepulse,nsamplepulse);

    // Symmetrize just in case
    //Eigen::MatrixXd block_sym = block_allocated.template triangularView<Eigen::Lower>();
    //block_sym = block_sym + block_sym.transpose().triangularView<Eigen::StrictlyUpper>();

    //std::cout << "block_allocated: " << block_allocated << std::endl;

    // Try LLT to check PD
    //Eigen::LLT<Eigen::MatrixXd> llt(block_allocated);
    //if (llt.info() == Eigen::NumericalIssue) {
    //    std::cout << "WARNING: pulse block not PD" << std::endl;
    //}

    _invcov.block(firstsamplet,firstsamplet,nsamplepulse,nsamplepulse).triangularView<Eigen::Lower>() += ampsq*block_allocated;
  }

  // std::cout << " updateCov " << " here "  << std::endl;
  //std::cout << " invcov after adding pulse covariance= " << std::endl << _invcov << std::endl;


  _invcov.triangularView<Eigen::Upper>() =
    _invcov.transpose().triangularView<Eigen::Upper>();



   // Quick PD check via LLT
   Eigen::LLT<SampleMatrix> llt(_invcov);
   if (llt.info() == Eigen::NumericalIssue) {
      std::cout << "WARNING: _invcov not PD" << std::endl;
   }

   _covdecomp.compute(_invcov);
    //  std::cout << " updateCov " << " done "  << std::endl;
  
  bool status = true;
  return status;
  
}

double PulseChiSqSNNLS::ComputeChiSq() {

    SampleVector model = _pulsemat*_ampvec;

    std::cout << "model pre-time: " << model << std::endl;
    // add timing shifts
    for (unsigned int ipulse=0; ipulse<_bxs.rows(); ++ipulse) {
        if (_timeActive(ipulse)) {
            std::cout << "ipulse (time Active here): " << ipulse << ", _time(ipulse): " << _time(ipulse) << ", _pulsemat_t.col(ipulse): " << _pulsemat_t.col(ipulse) << ", _ampvec(ipulse): " << _ampvec(ipulse) << std::endl;
            model -= _time(ipulse) * _pulsemat_t.col(ipulse) * _ampvec(ipulse);
        }
    }
    std::cout << "model post-time: " << model << std::endl;

    //std::cout << "_invcov at chi2 / residuals step: " << std::endl << _invcov << std::endl;

    //SampleMatrix L_debug = _covdecomp.matrixL();
    //std::cout << "_covdecomp.matrixL at chi2 / residuals step: " << std::endl << L_debug << std::endl;

    _normResVec = _covdecomp.matrixL().solve(model - _sampvec);
    _absResVec = model - _sampvec;
    //std::cout << "_normResVec: " << _normResVec << std::endl;
    //std::cout << "_absResVec: " << _absResVec << std::endl;
    return _normResVec.squaredNorm();
}

double PulseChiSqSNNLS::ComputeApproxUncertainty(unsigned int ipulse) {
  //compute approximate uncertainties
  //(using 1/second derivative since full Hessian is not meaningful in
  //presence of positive amplitude boundaries.)
  
  return 1./_covdecomp.matrixL().solve(_pulsemat.col(ipulse)).norm();
  
}

bool PulseChiSqSNNLS::NNLS() {
  
  //Fast NNLS (fnnls) algorithm as per http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.157.9203&rep=rep1&type=pdf
  
  const unsigned int npulse = _bxs.rows();
  
  SamplePulseMatrix invcovp = _covdecomp.matrixL().solve(_pulsemat);
  aTamat.triangularView<Eigen::Lower>() = invcovp.transpose()*invcovp;
  aTamat = aTamat.selfadjointView<Eigen::Lower>();
  aTbvec = invcovp.transpose()*_covdecomp.matrixL().solve(_sampvec);  
  
  
  PulseVector wvec(npulse);
  
  
  int iter = 0;
  while (true) {    
    //can only perform this step if solution is guaranteed viable
    if (iter>0 || _nP==0) {
      if ( _nP==npulse ) break;                  
      
      const unsigned int nActive = npulse - _nP;
      
      wvec.tail(nActive) = aTbvec.tail(nActive) - (aTamat.selfadjointView<Eigen::Lower>()*_ampvec).tail(nActive);       
      
      Index idxwmax;
      double wmax = wvec.tail(nActive).maxCoeff(&idxwmax);
      
      //convergence
      if (wmax<1e-11) break;

      //worst case protection
      if (iter >= 500) {
      	std::cout << "PulseChiSqSNNLS::NNLS()" << "\tMax Iterations reached at iter " << iter << std::endl;
        break;
      }
      
      //unconstrain parameter
      Index idxp = _nP + idxwmax;
      NNLSUnconstrainParameter(idxp);
      // std::cout << "adding index " << int(idxp) << " orig index " << int(_bxs.coeff(idxp)) << std::endl;
    }
    
    
    while (true) {
      // std::cout << "iter in, _nP = " << _nP << std::endl;
      
      // std::cout << " >>  iter = " << iter << std::endl;
      
      if (_nP==0) break;     
      
      ampvecpermtest = _ampvec;
      
      //solve for unconstrained parameters      
      //ampvecpermtest.head(_nP) = aTamat.topLeftCorner(_nP,_nP).ldlt().solve(aTbvec.head(_nP));     

      //need to have specialized function to call optimized versions
      // of matrix solver... this is truly amazing...
      eigen_solve_submatrix(aTamat, aTbvec, ampvecpermtest, _nP);
      
      //check solution
      bool positive = true;
      for (unsigned int i = 0; i < _nP; ++i)
        positive &= (ampvecpermtest(i) > 0);
      if (positive) {
        _ampvec.head(_nP) = ampvecpermtest.head(_nP);
        break;
      }
      
      //update parameter vector
      Index minratioidx=0;
      
      double minratio = std::numeric_limits<double>::max();
      for (unsigned int ipulse=0; ipulse<_nP; ++ipulse) {
        if (ampvecpermtest.coeff(ipulse)<=0.) {
	  const double c_ampvec = _ampvec.coeff(ipulse);
          const double ratio = c_ampvec/(c_ampvec - ampvecpermtest.coeff(ipulse));
          if (ratio<minratio) {
            minratio = ratio;
            minratioidx = ipulse;
          }
        }
      }
      
      _ampvec.head(_nP) += minratio*(ampvecpermtest.head(_nP) - _ampvec.head(_nP));
      
      //avoid numerical problems with later ==0. check
      _ampvec.coeffRef(minratioidx) = 0.;
      
      // std::cout << "removing index " << int(minratioidx) << " orig idx " << int(_bxs.coeff(minratioidx)) << std::endl;
      NNLSConstrainParameter(minratioidx);      
    }
    ++iter;
    
    
    //---- AM:: add this new check to stop
    if (iter > 1000) break;
    
  }
  
  // std::cout << "     -> _ampvec = " << std::endl << _ampvec << std::endl;
  
  return true;
  
  
}

void PulseChiSqSNNLS::NNLSUnconstrainParameter(Index idxp) {
  aTamat.col(_nP).swap(aTamat.col(idxp));
  aTamat.row(_nP).swap(aTamat.row(idxp));
  _pulsemat.col(_nP).swap(_pulsemat.col(idxp));
  _pulsemat_t.col(_nP).swap(_pulsemat_t.col(idxp));
  std::swap(aTbvec.coeffRef(_nP), aTbvec.coeffRef(idxp));
  std::swap(_ampvec.coeffRef(_nP), _ampvec.coeffRef(idxp));
  std::swap(_bxs.coeffRef(_nP), _bxs.coeffRef(idxp));
  ++_nP;
}

void PulseChiSqSNNLS::NNLSConstrainParameter(Index minratioidx) {
  aTamat.col(_nP - 1).swap(aTamat.col(minratioidx));
  aTamat.row(_nP - 1).swap(aTamat.row(minratioidx));
  _pulsemat.col(_nP - 1).swap(_pulsemat.col(minratioidx));
   _pulsemat_t.col(_nP - 1).swap(_pulsemat_t.col(minratioidx));
  std::swap(aTbvec.coeffRef(_nP - 1), aTbvec.coeffRef(minratioidx));
  std::swap(_ampvec.coeffRef(_nP - 1), _ampvec.coeffRef(minratioidx));
  std::swap(_bxs.coeffRef(_nP - 1), _bxs.coeffRef(minratioidx));
  --_nP;
}

bool PulseChiSqSNNLS::OnePulseMinimize() {
  //Fast NNLS (fnnls) algorithm as per http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.157.9203&rep=rep1&type=pdf

  //   const unsigned int npulse = 1;

  invcovp = _covdecomp.matrixL().solve(_pulsemat);
  //   aTamat = invcovp.transpose()*invcovp;
  //   aTbvec = invcovp.transpose()*_covdecomp.matrixL().solve(_sampvec);

  SingleMatrix aTamatval = invcovp.transpose() * invcovp;
  SingleVector aTbvecval = invcovp.transpose() * _covdecomp.matrixL().solve(_sampvec);
  _ampvec.coeffRef(0) = std::max(0., aTbvecval.coeff(0) / aTamatval.coeff(0));

  return true;
}
