#include "dg_external_state.h"
#include <vector>
#include "dg_user_defined.h"
#include <cmath>
#include "dg_param.h"
#include "dg_kx_ky.h"

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

void External_state_Gaussian_exact2(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){


	
	double kx{}, ky{};

	Get_kx_ky(kx, ky, t);

	double inter = exp( - pow((kx * (x - user::xx0) + 
			ky * (y - user::yy0) - dg_fun::C * t), 2) / (user::D * user::D) );

	q_ext[index[0]] = inter;

	q_ext[index[1]] = kx / dg_fun::C * inter;

	q_ext[index[2]] = ky / dg_fun::C * inter;
}

// test
void External_state_sin_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){

	double inter1 = std::sin(user::pi * x + t);
	double inter2 = std::sin(user::pi * y + t);

	q_ext[index[0]] = -1.0 / user::pi * inter1 - 1.0 / user::pi * inter2;

	q_ext[index[1]] = inter1;

	q_ext[index[2]] = inter2;
}
