#include <iostream>
#include "tools.h"
#include "math.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
	int len;
	VectorXd rmse(4);
	len = estimations.size();
	for(int i = 0; i < len; ++i)
	{
		VectorXd residual(len);
		residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	rmse << rmse/len;
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  float px = x_state[0];
  float py = x_state[1];
  float vx = x_state[2];
  float vy = x_state[3];
  MatrixXd jacobian(3, 4);

  if(fabs(pow(px, 2) + pow(py, 2)) < 0.0001)
  {
  	cout << "Error: division by zero" << endl;
  	return jacobian;
  }

  float a11, a12, a21, a22, a31, a32;
  a11 = px / sqrt(pow(px, 2) + pow(py, 2));
  a12 = py / sqrt(pow(px, 2) + pow(py, 2));
  a21 = (-py) / (pow(px, 2) + pow(py, 2));
  a22 = (px) / (pow(px, 2) + pow(py, 2));
  a31 = (py*((vx*py) - (vy*px))) / pow((pow(px, 2) + pow(py, 2)), 1.5);
  a32 = (px*((vy*px) - (vx*py))) / pow((pow(px, 2) + pow(py, 2)), 1.5);
  // a33 and a34 are same as a11 and a12

  jacobian << a11, a12, 0, 0,
  				a21, a22, 0, 0,
  				a31, a32, a11, a12;

  return jacobian;
}
