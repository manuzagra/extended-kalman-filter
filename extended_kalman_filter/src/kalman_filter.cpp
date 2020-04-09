#include "kalman_filter.h"

#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * predict the state
   */
	x_ = F_ * x_;
	P_ = (F_ * P_ * F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * update the state by using Kalman Filter equations
   */
	VectorXd y = z - (H_ * x_);
	MatrixXd S = (H_ * P_ * H_.transpose()) + R_;
    MatrixXd K =  P_ * H_.transpose() * S.inverse();

    // new state
    x_ = x_ + (K * y);
    P_ = (MatrixXd::Identity(x_.size(), x_.size()) - (K * H_)) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * update the state by using Extended Kalman Filter equations
   */
	// h(x)
	VectorXd h_x(3);
	h_x[0] = (double)sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
	h_x[1] = (double)atan2(x_[1], x_[0]);
	h_x[2] = (double)(x_[0]*x_[2] + x_[1]*x_[3]) / (sqrt(x_[0]*x_[0] + x_[1]*x_[1])); // / h_x[0]


	VectorXd y = z - h_x;
	while (y[1] >  M_PI) y[1] -= 2*M_PI;
	while (y[1] < -M_PI) y[1] += 2*M_PI;
	MatrixXd S = (H_ * P_ * H_.transpose()) + R_;
    MatrixXd K =  P_ * H_.transpose() * S.inverse();

    // new state
    x_ = x_ + (K * y);
    P_ = (MatrixXd::Identity(x_.size(), x_.size()) - (K * H_)) * P_;
}
