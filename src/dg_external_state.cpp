#include "dg_external_state.h"
#include <vector>
#include "dg_user_defined.h"
#include <cmath>
#include "dg_param.h"

// forward declaration -------------------------------------------------------------------------------------------------
void External_state_Gaussian_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_reflect_x(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index);

void External_state_reflect_y(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index);

void External_mirror_y_left(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);

void External_state_sin_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index);
//---------------------------------------------------------------------------------------------------------------------

/// @brief
/// Boundary conditions. User defined.
/// @param t current time.
/// @param x x coordinate (physical). 
/// @param y y coordinate (physical). 
/// @param q_ext exrernal solutions (outputs). 
/// @param index the index of three variables. 
void External_state_Gaussian_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){


	double inter = exp( - pow((user::kx * ( x - user::xx0) + 
			user::ky * (y - user::yy0) - dg_fun::C * t), 2) / (user::D * user::D) );

	q_ext[index[0]] = inter;

	q_ext[index[1]] = user::kx / dg_fun::C * inter;

	q_ext[index[2]] = user::ky / dg_fun::C * inter;
}

/// @brief
/// External state at bottom boundary comes from the mirror image wave.
/// Note that the domain should be [0, 1] * [0, 1]. 
/// @param t current time.
/// @param x x coordinate (physical). 
/// @param y y coordinate (physical). 
/// @param q_ext exrernal solutions (outputs). 
/// @param index the index of three variables. 
void External_mirror_y_left(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){

	static const double kx = - sqrt(2.0) / 2.0;
	static const double ky =   sqrt(2.0) / 2.0;
	static const double D = 0.2 / (2.0 * sqrt(log(2.0)));
	static const double x0 = 1.5;
	static const double y0 = 0.5;	

	double inter = exp( - pow((kx * (x - x0) + 
			ky * (y - y0) - dg_fun::C * t), 2) / (D * D) );

	q_ext[index[0]] = inter;

	q_ext[index[1]] = kx / dg_fun::C * inter;

	q_ext[index[2]] = ky / dg_fun::C * inter;

}



/// @brief
/// Reflect boundary conditions. 
void External_state_reflect_x(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index){

	q_ext[index[0]] = q_int[index[0]];
	
	q_ext[index[1]] = (2.0 * q_int[index[0]] + dg_fun::C * q_int[index[1]]) / dg_fun::C;

	q_ext[index[2]] = q_int[index[2]];

}

void External_state_reflect_y(std::vector<double>& q_int, std::vector<double>& q_ext, 
				std::vector<int>& index){

	q_ext[index[0]] = q_int[index[0]];

	q_ext[index[1]] = 0.0;

	q_ext[index[2]] = (2.0 * q_int[index[0]] + dg_fun::C * q_int[index[2]]) / dg_fun::C;


}


// test
void External_state_sin_exact(double t, double x, double y, std::vector<double>& q_ext, std::vector<int>& index){

	double inter1 = std::sin(user::pi * x + t);
	double inter2 = std::sin(user::pi * y + t);

	q_ext[index[0]] = -1.0 / user::pi * inter1 - 1.0 / user::pi * inter2;

	q_ext[index[1]] = inter1;

	q_ext[index[2]] = inter2;
}
