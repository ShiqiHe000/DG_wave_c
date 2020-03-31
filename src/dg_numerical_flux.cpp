#include <vector>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_external_state.h"
#include "dg_riemann_solver.h"
#include <algorithm>
#include "dg_nodal_2d_storage.h"
#include "dg_param.h"
#include "dg_affine_map.h"
#include "dg_numerical_flux.h"
#include <iostream>	//test

// forward declaration-----------------------------------------------------
void Numerical_flux_x(double t);
void Numerical_flux_y(double t);
//------------------------------------------------------------------------


/// @brief
/// Compute the numerical fluxes on the x direction interfaces for all the elements.  
/// @param t time.
void Numerical_flux_x(double t){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){
		
		int pordery = temp -> m;
		int size = dg_fun::num_of_equation * (pordery + 1);	// right now assume conforming interface

		temp -> nflux_l = std::vector<double>(size);
		temp -> nflux_r = std::vector<double>(size);

		// compute numerical flux on the south interface
		for(auto it_face = temp -> facen[0].begin(); it_face != temp -> facen[0].end(); ++it_face){

			// index for three current points (three equations).
			std::vector<int> index{0, pordery + 1, (pordery + 1) * 2};	

			if(it_face -> face_type == 'L'){	// local neighbour
				
				int n_key = it_face -> key;	// neighbour's key

				for(int s = 0; s <= pordery; ++s){

					//now conforming interface=========================================================
					// Riemann solver
					Riemann_solver_x(local::Hash_elem[n_key] -> solution_int_r, temp -> solution_int_l, 
							temp -> nflux_l, -1, index);
					//=================================================================================

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
					

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}


			}
		}

		if((temp -> ghost).size() > 0){	// if this element has remote neighbours, clear the ghost layer

			(temp -> ghost).clear();

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
				Riemann_solver_x(temp -> solution_int_r, solution_ext,
						temp -> nflux_r, 1, index);
				
				std::transform(index.begin(), index.end(), index.begin(), 
						[](int x){return (x + 1);});		// increment 1
			}
//if(mpi::rank == 1){
//
//	for(int i = 0; i < size; ++i){
//
//		std::cout<< "num " << i << " "<< (temp -> nflux_r)[i] << "\n";
//
//	}
//
//}

		}


		temp = temp -> next;
	}
}

/// @brief
/// Compute the numerical fluxes on the x direction interfaces for all the elements.  
/// @param t time.
void Numerical_flux_y(double t){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){
		
		int porderx = temp -> n;
		int size = dg_fun::num_of_equation * (porderx + 1);	// right now assume conforming interface

		temp -> nflux_l = std::vector<double>(size);
		temp -> nflux_r = std::vector<double>(size);

		// compute numerical flux on the west interface
		for(auto it_face = temp -> facen[2].begin(); it_face != temp -> facen[2].end(); ++it_face){

			// index for three current points (three equations).
			std::vector<int> index{0, porderx + 1, (porderx + 1) * 2};	

			if(it_face -> face_type == 'L'){	// local neighbour
				
				int n_key = it_face -> key;	// neighbour's key

				for(int s = 0; s <= porderx; ++s){

					// Riemann solver
					Riemann_solver_y(local::Hash_elem[n_key] -> solution_int_r, temp -> solution_int_l, 
							temp -> nflux_l, -1, index);


					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}	
			}
			else if(it_face -> face_type == 'B'){	// phsical boundary

				std::vector<double> solution_ext(size);

				double del_x = ((temp -> xcoords[1]) - (temp -> xcoords[0]));

				for(int s = 0; s <= porderx; ++s){
					// map to physical plane
					double x = Affine_mapping(nodal::gl_points[porderx][s], (temp -> xcoords[0]), del_x);

					// impose boundary conditions
					External_state_Gaussian_exact(t, x, (temp -> ycoords[0]), solution_ext, index);

					// Riemann solver
					Riemann_solver_y(solution_ext, temp -> solution_int_l, 
							temp -> nflux_l, -1, index);
					
					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}
			}
			else{	// mpi boundary

				int n_key = it_face -> key;

				for(int s = 0; s <= porderx; ++s){

					// Riemann solver
					Riemann_solver_y((temp -> ghost[n_key]), temp -> solution_int_l, 
							temp -> nflux_l, -1, index);
					

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}
			}
		}

		if((temp -> ghost).size() > 0){	// if this element has remote neighbours, clear the ghost layer

			(temp -> ghost).clear();

		}

		// compute numerical flux on the north interface (only elements face the physical boundary)
		auto it_face = temp -> facen[3].begin();
		if(it_face -> face_type == 'B'){

			std::vector<double> solution_ext(size);
	
			double del_x = ((temp -> xcoords[1]) - (temp -> xcoords[0]));
			
			// index for three current points (three equations).
			std::vector<int> index{0, porderx + 1, (porderx + 1) * 2};	
	
			for(int s = 0; s <= porderx; ++s){
				// map to physical plane
				double x = Affine_mapping(nodal::gl_points[porderx][s], (temp -> xcoords[0]), del_x);
	
				// impose boundary conditions
				External_state_Gaussian_exact(t, x, (temp -> ycoords[1]), solution_ext, index);
	
				// Riemann solver
				Riemann_solver_y(temp -> solution_int_r, solution_ext, 
						temp -> nflux_r, 1, index);
				
				std::transform(index.begin(), index.end(), index.begin(), 
						[](int x){return (x + 1);});		// increment 1
			}

		}


		temp = temp -> next;
	}
}
