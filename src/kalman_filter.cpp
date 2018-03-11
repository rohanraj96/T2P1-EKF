#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y_ = z - H_ * x_;
  MatrixXd S_ = H_ * P_ * (H_.transpose()) + R_;
  MatrixXd K_ = P_ * (H_.transpose()) * (S_.inverse());
  x_ = x_ + (K_ * y_);
  int rows = x_.rows();
  MatrixXd I = MatrixXd::Identity(rows, rows);

  P_ = (I - (K_ * H_)) * P_;
}

float phi_norm(float phi)
{
  float tan_theta = tan(phi);
  return tanh(tan_theta) * M_PI;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd x_transform = VectorXd(3);
  x_transform << 0, 0, 0;

  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  if(!(px == 0 && py == 0)){
    float rho = sqrt(pow(px, 2) + pow(py, 2));
    float phi = atan(py/px);
    float rho_dot = ((px * vx) + (py * vy)) / rho;
    x_transform[0] = rho;
    x_transform[1] = phi;
    x_transform[2] = rho_dot;
  }

  VectorXd y_ = z - x_transform;

  y_[1] = phi_norm(y_[1]);
  MatrixXd S_ = H_ * P_ * H_.transpose() + R_;
  MatrixXd K_ = P_ * H_.transpose() * S_.inverse();

  x_ = x_ + K_ * y_;
  int rows = x_.rows();
  MatrixXd I_ = MatrixXd::Identity(rows, rows);
  P_ = (I_ - (K_ * H_)) * P_;
}
