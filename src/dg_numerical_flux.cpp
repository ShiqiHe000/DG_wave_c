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
#include <functional>	// std::plus

// forward declaration-----------------------------------------------------
void Numerical_flux_x(double t);
void Numerical_flux_y(double t);
void Two_vectors_sum(std::vector<double>& a, std::vector<double>& b);
void Form_mortar_x(Unit* temp, int n_key);
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

				temp -> ghost[n_key] = std::vector<double> (size);	// store neighbour's solution in ghost

				Form_mortar_x(temp, n_key);

				// left element, form L2 projection matrix
//				if(it_face -> hlevel == (temp -> mortar.l_max) && )

				for(int s = 0; s <= pordery; ++s){

					//now functionally conforming interface==========================================
					// Riemann solver
					Riemann_solver_x(local::Hash_elem[n_key] -> solution_int_r, temp -> solution_int_l, 
							temp -> ghost[n_key], -1, index);
					//=================================================================================

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}	

				// add up the numerical flux 
				Two_vectors_sum(temp -> ghost[n_key], temp -> nflux_l);

			}
			else if(it_face -> face_type == 'B'){	// phsical boundary

				std::vector<double> solution_ext(size);

				double del_y = ((temp -> ycoords[1]) - (temp -> ycoords[0]));

				for(int s = 0; s <= pordery; ++s){
					// map to physical plane
					double y = Affine_mapping(nodal::gl_points[pordery][s], (temp -> ycoords[0]), del_y);

					// impose boundary conditions (wave) ------------------------------------------------
//					External_state_Gaussian_exact(t, temp -> xcoords[0], y, solution_ext, index);
					// ----------------------------------------------------------------------------------


					// test -----------------------------------------------------------------------------
					External_state_sin_exact(t, temp -> xcoords[0], y, solution_ext, index);
					// ----------------------------------------------------------------------------------

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

					// assuming functionally conforming------------------------------------------------
					// Riemann solver, 
					// use ghost layer to store the nflux_l, so that we can send the corresponding one
					// to its neighbour. 
					Riemann_solver_x((temp -> ghost[n_key]), temp -> solution_int_l, 
							temp -> ghost[n_key], -1, index);
					//---------------------------------------------------------------------------------

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}

				// add up the numerical flux 
				Two_vectors_sum(temp -> ghost[n_key], temp -> nflux_l);


			}
		}

//		if((temp -> ghost).size() > 0){	// if this element has remote neighbours, clear the ghost layer
//
//			(temp -> ghost).clear();
//
//		}

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

		}


		temp = temp -> next;
	}
}


/// @brief
/// Form mortar structure for x direciton.
/// @param temp pointer to the current element.
/// @param n_key neighbour's key.
void Form_mortar_x(Unit* temp, int n_key){

	temp -> mortar.n_max = std::max(local::Hash_elem[n_key] -> n, temp -> n);	// maximum poly order

	temp -> mortar.l_max = std::max(local::Hash_elem[n_key] -> index[2], temp -> index[2]);	// maximum level

	// left element coordinate mapping
	if((local::Hash_elem[n_key] -> index[2]) == (temp -> mortar.l_max)){	// left element is the smallest

		// smallest element's coordinate does not need to scale
		temp -> mortar.a_l = 0.0;
		temp -> mortar.b_l = 1.0;	
	}
	else{	// right element is smaller

		double p = - ((temp -> ref_x[0] + (temp -> ref_x[1]))) / 2.0;

		double q = 2.0 / ((temp -> ref_x[1] - temp -> ref_x[0]));

		double s_d = ((local::Hash_elem[n_key] -> ref_x[0]) + p) * q;
		double s_u = ((local::Hash_elem[n_key] -> ref_x[1]) + p) * q;

		temp -> mortar.a_l = (s_u + s_d) / 2.0;
		temp -> mortar.b_l = s_u - (temp -> mortar.a_l);
	}

	// right element coordinate mapping
	if((temp -> index[2]) == (temp -> mortar.l_max)){	// right element is the smallest

		temp -> mortar.a_r = 0.0;
		temp -> mortar.b_r = 1.0;	
		
	}
	else{	// left element is smaller

		double p = - ((local::Hash_elem[n_key] -> ref_x[0]) + (local::Hash_elem[n_key] -> ref_x[1])) / 2.0;

		double q = 2.0 / (local::Hash_elem[n_key] -> ref_x[1] - local::Hash_elem[n_key] -> ref_x[0]);

		double s_d = ((temp -> ref_x[0]) + p) * q;
		double s_u = ((temp -> ref_x[1]) + p) * q;

		temp -> mortar.a_r = (s_u + s_d) / 2.0;
		temp -> mortar.b_r = s_u - (temp -> mortar.a_r);

	}


	// allocate the space for psi
	int size = (temp -> mortar.n_max + 1) * dg_fun::num_of_equation; 	// 1d array for mortar
	temp -> mortar.psi_l = std::vector<double>(size);
	temp -> mortar.psi_r = std::vector<double>(size);
}


/// @brief
/// Compute the sum of two vectors a and b, and store the result in vector b. 
/// @note The two vectors should have same size.
/// @param a vector a.
/// @param b vector b.
void Two_vectors_sum(std::vector<double>& a, std::vector<double>& b){

	std::transform(b.begin(), b.end(), a.begin(), b.begin(), std::plus<double> ());

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

		// compute numerical flux on the east interface (only elements face the physical boundary)
		auto it_face = temp -> facen[3].begin();
		if(it_face -> face_type == 'B'){

			std::vector<double> solution_ext(size);
	
			double del_x = ((temp -> xcoords[1]) - (temp -> xcoords[0]));
			
			// index for three current points (three equations).
			std::vector<int> index{0, porderx + 1, (porderx + 1) * 2};	
	
			for(int s = 0; s <= porderx; ++s){
				// map to physical plane
				double x = Affine_mapping(nodal::gl_points[porderx][s], (temp -> xcoords[0]), del_x);
	
				// impose boundary conditions --------------------------------------------------------
//				External_state_Gaussian_exact(t, x, (temp -> ycoords[1]), solution_ext, index);
				//------------------------------------------------------------------------------------

				// test -------------------------------------------------------------------------------
				External_state_sin_exact(t, x, (temp -> ycoords[1]), solution_ext, index);
				//------------------------------------------------------------------------------------
	
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
