#include "dg_user_defined.h"
#include "dg_affine_map.h"
#include <cmath>
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"
#include "dg_single_index.h"
#include <unordered_map>
#include <vector>

// Global variables
static const double kx = sqrt(2.0) / 2.0; 
static const double ky = sqrt(2.0) / 2.0; 
static const double D = 0.2 / (2.0 * sqrt(log(2.0)));
static const double xx0 = 0.0;
static const double yy0 = 0.0; 


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

			double inter = exp( - pow((kx * (x - xx0) + 
						ky * (y - yy0) - dg_fun::C * t), 2) / (D * D) );
		
			int nodei = Get_single_index(i, j, m + 1);

			e[0][nodei] = inter;
			e[1][nodei]= kx / dg_fun::C *inter;
			e[2][nodei]= ky / dg_fun::C *inter;
		}

	}
	

}
