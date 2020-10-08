#ifndef DG_EXTERNAL_STATE_H
#define DG_EXTERNAL_STATE_H

#include <vector>
#include <cmath>

// global vars ---------------------------------------------------
namespace mirror{
	static const double kx = - sqrt(2.0) / 2.0;
	static const double ky =   sqrt(2.0) / 2.0;
	static const double D = 0.2 / (2.0 * sqrt(log(2.0)));
	static const double x0 = 1.5;
	static const double y0 = 0.5;	
}
//----------------------------------------------------------------

void External_state_Gaussian_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_sin_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_reflect_x(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index);

void External_mirror_right_boundary(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_reflect_y(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index);
#endif
