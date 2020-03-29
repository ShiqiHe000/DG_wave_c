#include <vector>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_external_state.h"
#include "dg_riemann_solver.h"
#include <algorithm>
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_affine_map.h"
#include <cassert>

// forward declaration-----------------------------------------------------
void Unary_minus(std::vector<double>& a, std::vector<double>& b);

void Numerical_flux_x(double t);
//------------------------------------------------------------------------


/// @brief
/// Compute the numerical fluxes on the x direction interfaces for all the elements.  
/// @param t time.
void Numerical_flux_x(double t){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){
		
		int pordery = temp -> m;
		int size = dg_fun::num_of_equation * pordery;	// right now assume conforming interface

		temp -> nflux_l = std::vector<double>(size);
		temp -> nflux_r = std::vector<double>(size);

		// compute numerical flux on the south interface
		for(auto it_face = temp -> facen[0].begin(); it_face != temp -> facen[0].end(); ++it_face){

			// index for three current points (three equations).
			std::vector<int> index{0, pordery + 1, (pordery + 1) * 2};	

			if(it_face -> face_type == 'L'){	// local neighbour
				
				int n_key = it_face -> key;	// neighbour's key

				for(int s = 0; s <= pordery; ++s){

					// Riemann solver
					Riemann_solver_x(local::Hash_elem[n_key] -> solution_int_r, temp -> solution_int_l, 
							temp -> nflux_l, -1, index);

					// numerical flux on the right interface (assume same porder)
//					Unary_minus(temp -> nflux_l, temp -> nflux_r);

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}	
			}
			else if(it_face -> face_type == 'B'){	// phsical boundary

				std::vector<double> solution_ext(size);

				double del_y = ((temp -> ycoords[1]) - (temp -> ycoords[0]));

				for(int s = 0; s <= pordery; ++s){
					// map to physical plane
					double y = Affine_mapping(nodal::gl_points[pordery][s], (temp -> ycoords[0]), del_y);

					// impose boundary conditions
					External_state_Gaussian_exact(t, temp -> xcoords[0], y, solution_ext, index);

					// Riemann solver
					Riemann_solver_x(solution_ext, temp -> solution_int_l, 
							temp -> nflux_l, -1, index);
					
					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}
			}
			else{	// mpi boundary

				int n_key = it_face -> key;

				for(int s = 0; s <= pordery; ++s){

					// Riemann solver
					Riemann_solver_x((temp -> ghost[n_key]), temp -> solution_int_l, 
							temp -> nflux_l, -1, index);
					
					// numerical flux on the right interface
//					Unary_minus(temp -> nflux_l, temp -> nflux_r);

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}
			}
		}

		// compute numerical flux on the north interface (only elements face the physical boundary)
		auto it_face = temp -> facen[1].begin();
		if(it_face -> face_type == 'B'){

			std::vector<double> solution_ext(size);
	
			double del_y = ((temp -> ycoords[1]) - (temp -> ycoords[0]));
			
			// index for three current points (three equations).
			std::vector<int> index{0, pordery + 1, (pordery + 1) * 2};	
	
			for(int s = 0; s <= pordery; ++s){
				// map to physical plane
				double y = Affine_mapping(nodal::gl_points[pordery][s], (temp -> ycoords[0]), del_y);
	
				// impose boundary conditions
				External_state_Gaussian_exact(t, temp -> xcoords[1], y, solution_ext, index);
	
				// Riemann solver
				Riemann_solver_x(solution_ext, temp -> solution_int_l, 
						temp -> nflux_l, 1, index);
				
				std::transform(index.begin(), index.end(), index.begin(), 
						[](int x){return (x + 1);});		// increment 1
			}

		}


		temp = temp -> next;
	}
}

/// @brief
/// Vector b = - vector a (same size)
/// @param a the vector to apply - operator.
/// @param b vector b gets the value.
void Unary_minus(std::vector<double>& a, std::vector<double>& b){

	assert(a.size() == b.size() && "The size of the two vector should be equal. ");

	std::transform(a.begin(), a.end(), b.begin(), [](double x){return -x;});
}
