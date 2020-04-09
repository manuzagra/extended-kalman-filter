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

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // measurement matrix - laser
  H_laser_ << 1., 0., 0., 0.,
              0., 1., 0., 0.;

  // measurement matrix - radar
  // it is non linear, so the Jacobian will be used instead

  // initilize ekf
  ekf_.x_ = VectorXd(4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);

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

    // first measurement
    cout << "EKF: " << endl;

	// save the timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // Initialize state.
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        ekf_.x_ << measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]),
        		measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]),
				0.,
				0.;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        ekf_.x_ << measurement_pack.raw_measurements_[0],
        		   measurement_pack.raw_measurements_[1],
				   0.,
				   0.;
    }
    // Initialize state covariance
    ekf_.P_ << 1., 0., 0.,    0.,
    		   0., 1., 0.,    0.,
			   0., 0., 1000., 0.,
			   0., 0., 0.,    1000.;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


  // increment in time
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  /**
   * Prediction
   */

  // state transition matrix F
  ekf_.F_ << 1., 0., dt, 0.,
		     0., 1., 0., dt,
			 0., 0., 1., 0.,
			 0., 0., 0., 1.;

  // process noise covariance matrix Q
  // Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  double noise_ax = 9.0;
  double noise_ay = 9.0;

  double dt2 = dt*dt;
  double dt3 = dt2*dt;
  double dt4 = dt3*dt;

  ekf_.Q_ <<  dt4/4*noise_ax, 0,              dt3/2*noise_ax, 0,
              0,              dt4/4*noise_ay, 0,              dt3/2*noise_ay,
              dt3/2*noise_ax, 0,              dt2*noise_ax,   0,
              0,              dt3/2*noise_ay, 0,              dt2*noise_ay;


  ekf_.Predict();

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  // measurement matrix
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  // measurement covariance matrix
	  ekf_.R_ =  R_radar_;
	  // update
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
	  // measurement matrix
	  ekf_.H_ = H_laser_;
	  // measurement covariance matrix
	  ekf_.R_ =  R_laser_;
	  // update
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
