#include "dg_user_defined.h"
#include "dg_affine_map.h"
#include <cmath>
#include "dg_param.h"
#include "dg_poly_level_and_order.h"
#include "dg_nodal_2d_storage.h"
#include "dg_single_index.h"

// Global variables
const double kx = sqrt(2.0) / 2.0; 
const double ky = sqrt(2.0) / 2.0; 
const double D = 0.2 / (2.0 * sqrt(log(2.0)));
const double x0 = 0.0;
const double y0 = 0.0; 


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
				double del_x, double del_y, double* e, double t){

	// get poly level
	int level_x = Poly_order_to_level(grid::nmin, n);
	int level_y = Poly_order_to_level(grid::nmin, m);

	for(int j = 0; j <= m; ++j ){

		double y = Affine_mapping(nodal::gl_p[level_y][j], y_d, del_y);
		
		for(int i = 0; i <= n; ++i){

			double x = Affine_mapping(nodal::gl_p[level_x][i], x_l, del_x);

			double inter = exp( - pow((kx * (x - x0) + 
						ky * (y - y0) - dg_fun::C * t), 2) / (D * D) );
		
			int nodei[dg_fun::num_of_equation]{};	
			for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
				nodei[equ] = Get_single_index_3d(i, j, equ, grid::nmin + 1, grid::nmin + 1);

			}		
			
			e[node[0]] = inter;
			e[node[1]] = kx / dg_fun::C *inter;
			e[node[2]] = ky / dg_fun::C *inter;
		}

	}
	

}
