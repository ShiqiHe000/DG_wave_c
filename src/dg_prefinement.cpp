#include "dg_prefinement.h"
#include <vector>
#include "dg_unit.h"
#include "dg_param.h"
#include "dg_interpolate_to_new_points.h"
#include "dg_nodal_2d_storage.h"
#include <unordered_map>
#include "dg_local_storage.h"
#include "dg_status_table.h"
#include "dg_cantor_pairing.h"	
#include <iostream>	// test

// forward declaration----------------------------------------
void p_refinement_apply(Unit* temp);

void p_refinement(int kt);

void Update_neighbours_facen(Unit* temp);
//------------------------------------------------------------

/// @brief
/// Apply p-refinement to the elements
void p_refinement(int kt){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		// test==============================================================================
//		int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
//
//		if(key == 16 && kt == 0){
//
//			p_refinement_apply(temp);
//
//		}
		//====================================================================================

		if(temp -> prefine && (temp -> n < grid::nmax)){	// if flag refine and not exceed the maximum poly order

			p_refinement_apply(temp);

			temp -> prefine = false;	// turn off the switch after refine the element

		}

		temp = temp -> next;
	}

}

/// @brief
/// Increase polynomial order with 2 if does not exceed the maximum polynomial order. 
/// Here we increase polynomial order in x and y direction together.
/// @param temp Pointer to the current element.
void p_refinement_apply(Unit* temp){

	int new_order = temp -> n + 2;

	// impose p-refinement if new_order is inside the range
	if(new_order <= grid::nmax){

		// interpolate solution to new GL points
		// form interpolation matrix in x and y direction (assume refinement in x and y driection)-------------
		std::vector<double> T;	// only need one interpolation matrix (porderx = pordery)
		Polynomial_interpolate_matrix(nodal::gl_points[temp -> n], nodal::gl_points[new_order], T);
		// ----------------------------------------------------------------------------------------------------

		std::unordered_map<int, std::vector<double>> middle;	// intermidiate matrix

		// y direction
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
			
			int start_old{};
			int start_new{};

			middle[equ] = std::vector<double> ((new_order + 1) * (temp -> n + 1));
			
			for(int i = 0; i <= (temp -> n); ++i){
	
				// interval == 1 since we are in teh y direction
				// restriction: children inderit parent's polynomial order. 
				Interpolate_to_new_points(new_order + 1, temp -> m + 1, T,
						temp -> solution[equ], middle[equ], start_old, start_new, 1);

				start_old += (temp -> m + 1);
				start_new += (new_order + 1);
			}

		}

//for(auto& v : temp -> solution[1]){
//
//	std::cout<< v << "\n";
//}
//std::cout<< "\n";

//for(auto& v :  middle[1]){
//
//	std::cout<< v << "\n";
//}
//std::cout<< "\n";

		// clear the old solutions
		temp -> solution.clear();

		// x direction
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
			int start_old{};
			int start_new{};
			
			// allocate space
			temp -> solution[equ] = std::vector<double> ((new_order + 1) * (new_order + 1));

			for(int j = 0; j <= new_order; ++j){
	
				Interpolate_to_new_points(new_order + 1,  temp -> n + 1, T,
						middle[equ], temp -> solution[equ], start_old, start_new, new_order + 1);
				
				++start_old;
				++start_new;
			}

		}
	
//for(auto& v : temp -> solution[1]){
//
//	std::cout<< v << "\n";
//}
//std::cout<< "\n";
	
		// update the polynomial order
		temp -> n = new_order;
		temp -> m = new_order;

		// updates local neighbours face info
		Update_neighbours_facen(temp);

	}	
}

/// @brief
/// Since the current element's polynomial order changes, we need to update it's neighbours records.
/// @param temp Pointer to the current element.
void Update_neighbours_facen(Unit* temp){

	int my_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

	// loop through four faces
	for(int i = 0; i < 4; ++i){

		auto it_face = temp -> facen[i].begin();
		
		for(; it_face != temp -> facen[i].end(); ++it_face){
	
			int opp_dir = Opposite_dir(i);

			// only change local neighbours
			if(it_face -> face_type == 'L'){

				// go to neighbour and updates the orders
				for(auto it_n = local::Hash_elem[it_face -> key] -> facen[opp_dir].begin();
					it_n != local::Hash_elem[it_face -> key] -> facen[opp_dir].end(); ++it_n ){

					if(it_n -> key == my_key){	// find it

						it_n -> porderx = temp -> n;
						it_n -> pordery = temp -> m;
					}
				}
			}
			
		}

	}

}
