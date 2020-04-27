#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_param.h"
#include "dg_affine_map.h"
#include "dg_nodal_2d_storage.h"
#include <cmath>	// exp
#include "dg_single_index.h"
#include "dg_user_defined.h"
#include <iostream>	// test

/// @brief
/// Initialization all local elements based on the initial conditions.
void DG_init(){
	
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
			
				int num_p = Get_single_index(i, j, grid::nmin + 1);

				// wave --------------------------------------------------------------------------------	
				double inter = exp( - std::pow((user::kx * (x - user::xx0) + 
							user::ky * (y - user::yy0)), 2) / std::pow(user::D, 2));
      			
				temp -> solution[0][num_p] = inter;
				temp -> solution[1][num_p] = user::kx / dg_fun::C * inter;
				temp -> solution[2][num_p] = user::ky / dg_fun::C * inter;
				// -------------------------------------------------------------------------------------
				

				// test case sin(pi * x) + sin(pi * y) --------------------------------------------------
//				double inter1 = std::sin(user::pi * x);
//				double inter2 = std::sin(user::pi * y);
//
//				temp -> solution[0][num_p] = - 1.0 / user::pi * inter1 - 1.0 / user::pi * inter2;
//				temp -> solution[1][num_p] = inter1;
//				temp -> solution[2][num_p] = inter2;
//				
				//----------------------------------------------------------------------------------------

			}
		}
	
		temp = temp -> next;		

	}


}
