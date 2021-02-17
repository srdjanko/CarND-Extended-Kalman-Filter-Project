#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
   using namespace std;

  // check the validity of the following inputs:
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size())
  {
    throw "Estimation vector size must be equal to ground truth vector size.";
  }
  const int size = estimations.size();

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  VectorXd sum(4);

  for (int i=0; i < size; ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // calculate mean
  rmse = rmse / size;

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */

  // Jakobian is calculated as a part of 'void FusionEKF::ProcessMeasurement(...)'
  // since calculations are also used for measurement function (h).
}
