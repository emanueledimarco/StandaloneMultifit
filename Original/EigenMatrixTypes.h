#include <Eigen/Dense>

typedef Eigen::Matrix<double,16,1> SampleVector;
typedef Eigen::Matrix<double,20,1> FullSampleVector;
typedef Eigen::Matrix<double,Eigen::Dynamic,1,0,16,1> PulseVector;
typedef Eigen::Matrix<char,Eigen::Dynamic,1,0,16,1> BXVector;
typedef Eigen::Matrix<double,16,16> SampleMatrix;
typedef Eigen::Matrix<double,20,20> FullSampleMatrix;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,16,16> PulseMatrix;
typedef Eigen::Matrix<double,16,Eigen::Dynamic,0,16,16> SamplePulseMatrix;
typedef Eigen::LLT<SampleMatrix> SampleDecompLLT;

