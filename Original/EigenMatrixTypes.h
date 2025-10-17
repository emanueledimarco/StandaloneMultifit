//#include <Eigen/Dense>
#include </usr/include/eigen3/Eigen/Dense>

int SampleVectorSize = 10;
int FullSampleVectorSize = 19;
int PulseVectorSize = 12;

typedef Eigen::Matrix<double,SampleVectorSize,1> SampleVector;
typedef Eigen::Matrix<double,FullSampleVectorSize,1> FullSampleVector;
typedef Eigen::Matrix<double,Eigen::Dynamic,1,0,PulseVectorSize,1> PulseVector;
typedef Eigen::Matrix<char,Eigen::Dynamic,1,0,PulseVectorSize,1> BXVector;
typedef Eigen::Matrix<double,PulseVectorSize,PulseVectorSize> SampleMatrix;
typedef Eigen::Matrix<double,FullSampleVectorSize,FullSampleVectorSize> FullSampleMatrix;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,PulseVectorSize,PulseVectorSize> PulseMatrix;
typedef Eigen::Matrix<double,PulseVectorSize,Eigen::Dynamic,0,PulseVectorSize,PulseVectorSize> SamplePulseMatrix;
typedef Eigen::LLT<SampleMatrix> SampleDecompLLT;
