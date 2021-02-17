#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // state vector
  VectorXd x(4);
  x << 0,0,0,0;

  // transition matrix (without dt)
  MatrixXd F(4,4);
  F << 1,  0,  0,  0,
       0,  1,  0,  0,
       0,  0,  1,  0,
       0,  0,  0,  1;

  // state covariance matrix P
  MatrixXd P(4, 4);
  P << 10, 0, 0, 0,
       0, 10, 0, 0,
       0, 0, 100, 0,
       0, 0, 0, 100;

  // measurement matrix for laser measurements
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Jakobian matrix for radar measurements
  Hj_ = MatrixXd(3, 4);

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  // measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // process covariance matrix
  MatrixXd Q(4,4);

  // set the acceleration noise components
  noise_ax = 9;
  noise_ay = 9;

  ekf_.Init(x, P, F, H_laser_, R_laser_, Q);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      auto rho = measurement_pack.raw_measurements_[0];
      auto theta = measurement_pack.raw_measurements_[1];
      auto rho_dot = measurement_pack.raw_measurements_[2];
      // If we assume that theta_dot = 0, then we are still left
      // with a part of transformation which depends on {theta, rho_dot}
      // only.
      ekf_.x_ << rho * cos(theta),
                 rho * sin(theta),
                 rho_dot * cos(theta),
                 rho_dot * sin(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                measurement_pack.raw_measurements_[1],
                0,
                0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_2 * dt_2;

  ekf_.Q_ << dt_4 * noise_ax / 4, 0, dt_3 * noise_ax / 2, 0,
              0,  dt_4 * noise_ay / 4, 0, dt_3 * noise_ay / 2,
              dt_3 * noise_ax / 2, 0, dt_2 * noise_ax, 0,
              0, dt_3 * noise_ay / 2, 0, dt_2 * noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // TODO: Radar updates
    ekf_.R_ = R_radar_;

    VectorXd z(3);
    z << measurement_pack.raw_measurements_[0],
         measurement_pack.raw_measurements_[1],
         measurement_pack.raw_measurements_[2];

    auto x = ekf_.x_[0]; auto y = ekf_.x_[1];
    auto x_dt = ekf_.x_[2]; auto y_dt = ekf_.x_[3];
    auto m1 = x*x + y*y;

    VectorXd h(3);
    MatrixXd Hj(3, 4);

    // check division by zero
    if (fabs(m1) < 0.0001) {
      cout << "CalculateJacobian () - Error - Division by Zero" << endl;
      return;
    }
    else
    {
      // I calculate here the angle (phi) difference between measurement and estimation
      // by using the quotient of their complex counterparts:
      // (cos_phi + I * sin_phi) / (x + I * y).
      // The angle difference can be obtained by applying the atan2 to the components
      // of the quotient.

      const auto cos_phi = cos(z[1]);
      const auto sin_phi = sin(z[1]);
      auto angle_diff = atan2(-y * cos_phi + x * sin_phi, x * cos_phi + y * sin_phi);
      z[1] = angle_diff;

      auto m2 = sqrt(m1);
      auto m3 = x * x_dt + y * y_dt;
      auto m4 = m1 * m2;

      h << m2,
            0, // atan2(y, x),    // The difference is set in z[1]
            m3 / m2;

      Hj << x / m2, y / m2, 0, 0,
            -y / m1, x / m1, 0, 0,
            y * (x_dt * y - x * y_dt) / m4, x * (x * y_dt - x_dt * y) / m4, x / m2, y / m2;
    }

    ekf_.UpdateEKF(z, h, Hj);
  }
  else {

    // TODO: Laser updates
    ekf_.R_ = R_laser_;

    VectorXd z(2);
    z << measurement_pack.raw_measurements_[0],
         measurement_pack.raw_measurements_[1];

    ekf_.Update(z);
  }

  // print the output
//  cout << "x_ = " << ekf_.x_ << endl;
//  cout << "P_ = " << ekf_.P_ << endl;
}
