#include "dg_external_state.h"
#include <vector>
#include "dg_user_defined.h"
#include <cmath>
#include "dg_param.h"

/// @brief
/// Boundary conditions. User defined.
/// @param t current time.
/// @param x x coordinate (physical). 
/// @param y y coordinate (physical). 
/// @param 
void External_state_Gaussian_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){


	double inter = exp( - pow((user::kx * (x - user::xx0) + 
			user::ky * (y - user::yy0) - dg_fun::C * t), 2) / (user::D * user::D) );

	q_ext[index[0]] = inter;

	q_ext[index[1]] = user::kx / dg_fun::C * inter;

	q_ext[index[2]] = user::ky / dg_fun::C * inter;
}
