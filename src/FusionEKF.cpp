#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  noise_ax = 9.0;
  noise_ay = 9.0;

  R_laser_ << 0.0225, 0,
        0, 0.0225;

  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

  ekf_.Q_ = MatrixXd(4, 4);

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  if (!is_initialized_) {

    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      const float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      ekf_.x_ << rho * cos(phi), rho * sin(phi), rho_dot * cos(phi), rho_dot * sin(phi);
      previous_timestamp_ = measurement_pack.timestamp_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      previous_timestamp_ = measurement_pack.timestamp_;
    }

    is_initialized_ = true;
    return;
  }

  // cout << "CALCULATE DEL_T" << endl;
  float del_t = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;

  if(del_t > 0.001){

  previous_timestamp_ = measurement_pack.timestamp_;
  // cout << "TIME: " << del_t << endl;

  ekf_.F_ << 1, 0, del_t, 0,
              0, 1, 0, del_t,
              0, 0, 1, 0,
              0, 0, 0, 1;

  float dt_4 = pow(del_t, 4);
  float dt_3 = pow(del_t, 3);
  float dt_2 = pow(del_t, 2);

  ekf_.Q_ << (dt_4 * noise_ax) / 4, 0, (dt_3 * noise_ax) / 2, 0,
              0, (dt_4 * noise_ay) / 4, 0, (dt_3 * noise_ay) / 2,
              (dt_3 * noise_ax) / 2, 0, (dt_2 * noise_ax), 0,
              0, (dt_3 * noise_ay) / 2, 0, (dt_2 * noise_ay);

  ekf_.Predict();


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    cout << "RADAR UPDATE" << endl;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  else {

    cout << "LASER UPDATE" << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout << "END OF FUSIONEKF.CPP" << endl;
}
}