#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_param.h"
#include "dg_affine_map.h"
#include "dg_nodal_2d_storage.h"
#include <cmath>	// exp
#include "dg_single_index.h"

/// @brief
/// Initialization all local elements based on the initial conditions.
void DG_init(){
	
	// coefficients in the I.Cs
	const double kx = sqrt(2.0) / 2.0;
	const double ky = sqrt(2.0) / 2.0;
	const double D = 0.2 / (2.0 * sqrt(log(2.0)));
	const double x0 = 0.0;
	const double y0 = 0.0;


	Unit* temp = local::head;

	// traverse the linked list
	for(int k = 0; k < local::local_elem_num; ++k){

		// elemement size
		double del_x = (temp -> xcoords[1]) - (temp -> xcoords[0]);
		double del_y = (temp -> ycoords[1]) - (temp -> ycoords[0]);

		for(int j = 0; j <= grid::nmin; ++j){
			
			// map reference location to physical localtion
			double gl_p_y = nodal::gl_points[grid::nmin][j];
			double y = Affine_mapping(gl_p_y, temp -> ycoords[0], del_y);
			
			for(int i = 0; i <= grid::nmin; ++i){
	
				double gl_p_x = nodal::gl_points[grid::nmin][i];
				double x = Affine_mapping(gl_p_x, temp -> xcoords[0], del_x);
				
				double inter = exp( - std::pow((kx * (x - x0) + 
							ky * (y - y0)), 2) / std::pow(D, 2));
				
				int num_p = Get_single_index(i, j, grid::nmin + 1);

				temp -> solution[0][num_p] = inter;
				temp -> solution[1][num_p] = kx / dg_fun::C * inter;
				temp -> solution[2][num_p] = ky / dg_fun::C * inter;
			}
		}
	
		temp = temp -> next;		

	}


}
