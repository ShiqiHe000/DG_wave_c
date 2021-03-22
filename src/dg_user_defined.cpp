#include "dg_user_defined.h"
#include "dg_affine_map.h"
#include <cmath>
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"
#include "dg_single_index.h"
#include <vector>
#include <unordered_map>
#include "dg_kx_ky.h"

// Global variables
namespace user{
	const double kx = sqrt(2.0) / 2.0; 
	const double ky = - sqrt(2.0) / 2.0; 
//	const double kx = 0.0; 
//	const double ky = 1.0; 
	const double D = 0.2 / (2.0 * sqrt(log(2.0)));
	const double xx0 = 0.5;
	const double yy0 = 0.5; 

	const double pi = 4.0 * std::atan(1.0);	// testing

};


/// @brief
/// Exact solution  
/// @param n poly order in x driection
/// @param m poly order in y direction
/// @param gl_x current GL point in x direction
/// @param gl_y current GL point in y direction
/// @param x_l element left boundary in x direction (reference)
/// @param y_d element left boundary in y direction (reference)
/// @param del_x element size in x direction
/// @param del_y element size in y direction
/// @param e exact solution
/// @param t current time step
void Exact_solution_Gaussian(int n, int m, double x_l, double y_d,
				double del_x, double del_y, std::unordered_map<int, std::vector<double>>& e, double t){

	for(int j = 0; j <= m; ++j ){

		double y = Affine_mapping(nodal::gl_points[m][j], y_d, del_y);
		
		for(int i = 0; i <= n; ++i){

			double x = Affine_mapping(nodal::gl_points[n][i], x_l, del_x);

			double inter = exp( - pow((user::kx * (x - user::xx0) + 
						user::ky * (y - user::yy0) - dg_fun::C * t), 2) / (user::D * user::D) );
		
			int nodei = Get_single_index(i, j, m + 1);

			e[0][nodei] = inter;
			e[1][nodei]= user::kx / dg_fun::C * inter;
			e[2][nodei]= user::ky / dg_fun::C * inter;
		}

	}
	

}

void Exact_solution_Gaussian2(int n, int m, double x_l, double y_d,
				double del_x, double del_y, std::unordered_map<int, std::vector<double>>& e, double t){

	for(int j = 0; j <= m; ++j ){

		double y = Affine_mapping(nodal::gl_points[m][j], y_d, del_y);
		
		for(int i = 0; i <= n; ++i){

			double x = Affine_mapping(nodal::gl_points[n][i], x_l, del_x);

			double kx{}, ky{};

			Get_kx_ky(kx, ky, t);
		
			double inter = exp( - pow((kx * (x - user::xx0) + 
						ky * (y - user::yy0) 
						- dg_fun::C * t), 2) / (user::D * user::D) );
		
			int nodei = Get_single_index(i, j, m + 1);

			e[0][nodei] = inter;
			e[1][nodei]= kx / dg_fun::C * inter;
			e[2][nodei]= ky / dg_fun::C * inter;
		}

	}
	

}

// test case
void Exact_solution_sin(int n, int m, double x_l, double y_d,
				double del_x, double del_y, std::unordered_map<int, std::vector<double>>& e, double t){

	for(int j = 0; j <= m; ++j ){

		double y = Affine_mapping(nodal::gl_points[m][j], y_d, del_y);
		
		for(int i = 0; i <= n; ++i){

			double x = Affine_mapping(nodal::gl_points[n][i], x_l, del_x);

			double inter1 = std::sin(user::pi * x + t);
			double inter2 = std::sin(user::pi * y + t);
	
			int nodei = Get_single_index(i, j, m + 1);

			e[0][nodei] = - 1.0 / user::pi * inter1 - 1.0 / user::pi * inter2;
			e[1][nodei] = inter1;
			e[2][nodei] = inter2;
		}

	}
	

}
