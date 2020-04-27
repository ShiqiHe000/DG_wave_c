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
#include "dg_mortar_construction.h"
#include <iostream>	// test

// forward declaration-----------------------------------------------------
void Numerical_flux_x(double t);
void Numerical_flux_y(double t);
void Two_vectors_sum(std::vector<double>& a, std::vector<double>& b);	// useless
void Form_mortar_x(Unit* temp, const std::vector<Unit::Face>::iterator it_face);
void Form_mortar_y(Unit* temp, const std::vector<Unit::Face>::iterator it_face);
void Gen_a_and_b(double zd, double zu, double sd, double su, double& a, double& b);
//------------------------------------------------------------------------


	/// @brief
/// Compute the numerical fluxes on the x direction interfaces for all the elements.  
/// @param t time.
void Numerical_flux_x(double t){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){
		
		int pordery = temp -> m;
		int size = dg_fun::num_of_equation * (pordery + 1);	// right element interface

		temp -> nflux_l = std::vector<double>(size);
		temp -> nflux_r = std::vector<double>(size);

		// compute numerical flux on the south interface
		for(auto it_face = temp -> facen[0].begin(); it_face != temp -> facen[0].end(); ++it_face){


			if(it_face -> face_type == 'L'){	// local neighbour
				
				int n_key = it_face -> key;	// neighbour's key

				int size_n = (it_face -> pordery + 1) * dg_fun::num_of_equation; // left element interface

				temp -> ghost[n_key] = std::vector<double> (size_n);	// store neighbour's solution in ghost

				Form_mortar_x(temp, it_face); // allocate space on the mortar

				std::vector<double> T;	// interpolaiton matrix

				// left element, L2 projection
				L2_projection_to_mortar(temp -> mortar.n_max, it_face -> pordery,
								it_face -> hlevel, temp -> mortar.l_max, 
								temp -> mortar.a_l, temp -> mortar.b_l,
					 			local::Hash_elem[n_key] -> solution_int_r, 
								temp -> mortar.psi_l, T);

				// right element, L2 projection
				L2_projection_to_mortar(temp -> mortar.n_max, temp -> m,
								temp -> index[2], temp -> mortar.l_max, 
								temp -> mortar.a_r, temp -> mortar.b_r,
						 		temp -> solution_int_l, 
								temp -> mortar.psi_r, T);
				std::vector<int> index{0, (temp -> mortar.n_max + 1), (temp -> mortar.n_max + 1) * 2};	

				for(int s = 0; s <= temp -> mortar.n_max; ++s){

					//now functionally conforming interface==========================================
					// Riemann solver
					Riemann_solver_x(temp -> mortar.psi_l, temp -> mortar.psi_r, 
							temp -> mortar.nflux, -1, index);
					//=================================================================================

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}	

				// L2 project back to element, right element
				L2_projection_to_element(temp -> mortar.n_max, temp -> m, 
							temp -> index[2], temp -> mortar.l_max, 
							temp -> mortar.b_r,
			 				temp -> nflux_l, temp -> mortar.nflux, T);

				// L2 projection from mortar to left element	
				L2_projection_to_element(temp -> mortar.n_max, it_face -> pordery, 
							it_face -> hlevel, temp -> mortar.l_max, 
							temp -> mortar.b_l,
			 				temp -> ghost[n_key], temp -> mortar.nflux, T);


			}
			else if(it_face -> face_type == 'B'){	// phsical boundary

				std::vector<double> solution_ext(size);

				double del_y = ((temp -> ycoords[1]) - (temp -> ycoords[0]));

				std::vector<int> index{0, (pordery + 1), (pordery + 1) * 2};	

				for(int s = 0; s <= pordery; ++s){
					// map to physical plane
					double y = Affine_mapping(nodal::gl_points[pordery][s], (temp -> ycoords[0]), del_y);

					// impose boundary conditions (wave) ------------------------------------------------
					External_state_Gaussian_exact(t, temp -> xcoords[0], y, solution_ext, index);
					// ----------------------------------------------------------------------------------


					// test -----------------------------------------------------------------------------
					//External_state_sin_exact(t, temp -> xcoords[0], y, solution_ext, index);
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
				
				Form_mortar_x(temp, it_face); // allocate space on the mortar

				std::vector<double> T;	// only one side needs to be mapped

				// left element, L2 projection
				L2_projection_to_mortar(temp -> mortar.n_max, it_face -> pordery,
								it_face -> hlevel, temp -> mortar.l_max, 
								temp -> mortar.a_l, temp -> mortar.b_l,
						 		temp -> ghost[n_key], 
								temp -> mortar.psi_l, T);

				// right element, L2 projection
				L2_projection_to_mortar(temp -> mortar.n_max, temp -> m,
								temp -> index[2], temp -> mortar.l_max, 
								temp -> mortar.a_r, temp -> mortar.b_r,
						 		temp -> solution_int_l, 
								temp -> mortar.psi_r, T);
//if(mpi::rank == 1){
//
//	std::cout<< n_key << "\n";
//	std::cout.precision(17);
//	for(auto& h : temp -> mortar.psi_r){
//
//		std::cout<< std::fixed <<h << "\n";
//
//	}
//	std::cout<< "\n";
//
//}

				std::vector<int> index{0, (temp -> mortar.n_max + 1), (temp -> mortar.n_max + 1) * 2};	

				for(int s = 0; s <= temp -> mortar.n_max; ++s){

					// assuming functionally conforming------------------------------------------------
					// Riemann solver, 
					// use ghost layer to store the nflux_l, so that we can send the corresponding one
					// to its neighbour. 
					Riemann_solver_x(temp -> mortar.psi_l, temp -> mortar.psi_r, 
							temp -> mortar.nflux, -1, index);
					//---------------------------------------------------------------------------------
					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}

				// L2 project back to element, right element
				L2_projection_to_element(temp -> mortar.n_max, temp -> m, 
							temp -> index[2], temp -> mortar.l_max, 
							temp -> mortar.b_r,
			 				temp -> nflux_l, temp -> mortar.nflux, T);


//if(mpi::rank == 1){
//
//	std::cout<< n_key << "\n";
//	std::cout.precision(17);
//	for(auto& h : temp -> nflux_l){
//
//		std::cout<< std::fixed <<h << "\n";
//
//	}
//	std::cout<< "\n";
//
//}
				// L2 projection from mortar to left element	
				// store remote element's nunerical flux in ghost layer. But first clean up ghost layer.
				std::fill(temp -> ghost[n_key].begin(), temp -> ghost[n_key].end(), 0);

				L2_projection_to_element(temp -> mortar.n_max, it_face -> pordery, 
							it_face -> hlevel, temp -> mortar.l_max, 
							temp -> mortar.b_l,
			 				temp -> ghost[n_key], temp -> mortar.nflux, T);

//if(mpi::rank == 1){
//
//	std::cout<< n_key << "\n";
//	std::cout.precision(17);
//	for(auto& h : temp -> mortar.nflux){
//
//		std::cout<< std::fixed <<h << "\n";
//
//	}
//	std::cout<< "\n";
//
//}
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
	
				// impose boundary conditions (test) ------------------------------------------------
				//External_state_sin_exact(t, temp -> xcoords[1], y, solution_ext, index);
				// ----------------------------------------------------------------------------------

				// impose boundary conditions--------------------------------------------------------
				External_state_Gaussian_exact(t, temp -> xcoords[1], y, solution_ext, index);
				// ----------------------------------------------------------------------------------
	
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
/// @param temp pointer to the current element (right element).
/// @param it_face iterator points to the neighbour info in facen (left element).
void Form_mortar_x(Unit* temp, const std::vector<Unit::Face>::iterator it_face){

	// first clear the former info
	temp -> mortar = {}; 	// reset

	temp -> mortar.n_max = std::max(it_face -> pordery, temp -> m);	// maximum poly order

	temp -> mortar.l_max = std::max(it_face -> hlevel, temp -> index[2]);	// maximum level

	double zd, zu;	// mortar coord

	// left element coordinate mapping
	if((it_face -> hlevel) == (temp -> mortar.l_max)){	// left element is the smallest

		// smallest element's coordinate does not need to scale
		temp -> mortar.a_l = 0.0;
		temp -> mortar.b_l = 1.0;	

		zd = it_face -> ref_y[0];
		zu = it_face -> ref_y[1];

	}
	else{	// right element is smaller

		zd = temp -> ref_y[0];
		zu = temp -> ref_y[1];

		Gen_a_and_b(zd, zu, it_face -> ref_y[0], it_face -> ref_y[1], temp -> mortar.a_l, temp -> mortar.b_l);
	
	}

	// right element coordinate mapping
	if((temp -> index[2]) == (temp -> mortar.l_max)){	// right element is the smallest

		temp -> mortar.a_r = 0.0;
		temp -> mortar.b_r = 1.0;	
		
	}
	else{	// left element is smaller

		Gen_a_and_b(zd, zu, temp -> ref_y[0], temp -> ref_y[1], temp -> mortar.a_r, temp -> mortar.b_r);

	}


	// allocate the space for psi
	int size = (temp -> mortar.n_max + 1) * dg_fun::num_of_equation; 	// 1d array for mortar
	temp -> mortar.psi_l = std::vector<double>(size);
	temp -> mortar.psi_r = std::vector<double>(size);
	temp -> mortar.nflux = std::vector<double>(size);
}

/// @brief
/// Generate the offset a and scaling coefficient b.
/// @param zd lower boundary of the mortar element.
/// @param zu Upper boundary of the mortar element.
/// @param sd Upper boundary of the actual element.
/// @param su Upper boundary of the actual element.
/// @param a offset a.
/// @param b scaling ceofficient b. 
void Gen_a_and_b(double zd, double zu, double sd, double su, double& a, double& b){

	b = (sd - su) / (zd - zu);

	a = sd - zd * b;

}


/// @brief
/// Form mortar structure for y direciton.
/// @param temp pointer to the current element (right element).
/// @param it_face iterator points to the neighbour info in facen (left element).
void Form_mortar_y(Unit* temp, const std::vector<Unit::Face>::iterator it_face){

	// first clear the former info
	temp -> mortar = {}; 	// reset

	temp -> mortar.n_max = std::max(it_face -> porderx, temp -> n);	// maximum poly order

	temp -> mortar.l_max = std::max(it_face -> hlevel, temp -> index[2]);	// maximum level

	double zd, zu;
	
	// left element coordinate mapping
	if((it_face -> hlevel) == (temp -> mortar.l_max)){	// left element is the smallest

		// smallest element's coordinate does not need to scale
		temp -> mortar.a_l = 0.0;
		temp -> mortar.b_l = 1.0;	

		zd = it_face -> ref_x[0];
		zu = it_face -> ref_x[1];
	}
	else{	// right element is smaller
		zd = temp -> ref_x[0];
		zu = temp -> ref_x[1];

		Gen_a_and_b(zd, zu, it_face -> ref_x[0], it_face -> ref_x[1], temp -> mortar.a_l, temp -> mortar.b_l);
	}

	// right element coordinate mapping
	if((temp -> index[2]) == (temp -> mortar.l_max)){	// right element is the smallest

		temp -> mortar.a_r = 0.0;
		temp -> mortar.b_r = 1.0;	
		
	}
	else{	// left element is smaller

		Gen_a_and_b(zd, zu, temp -> ref_x[0], temp -> ref_x[1], temp -> mortar.a_r, temp -> mortar.b_r);
	}


	// allocate the space for psi
	int size = (temp -> mortar.n_max + 1) * dg_fun::num_of_equation; 	// 1d array for mortar
	temp -> mortar.psi_l = std::vector<double>(size);
	temp -> mortar.psi_r = std::vector<double>(size);
	temp -> mortar.nflux = std::vector<double>(size);
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

			if(it_face -> face_type == 'L'){	// local neighbour
				
				int n_key = it_face -> key;	// neighbour's key

				int size_n = (it_face -> pordery + 1) * dg_fun::num_of_equation; // left element interface

				temp -> ghost[n_key] = std::vector<double> (size_n);	// store neighbour's solution in ghost

				Form_mortar_y(temp, it_face); // allocate space on the mortar

				std::vector<double> T;

				// left element, L2 projection
				L2_projection_to_mortar(temp -> mortar.n_max, it_face -> porderx,
								it_face -> hlevel, temp -> mortar.l_max, 
								temp -> mortar.a_l, temp -> mortar.b_l,
						 		local::Hash_elem[n_key] -> solution_int_r, 
								temp -> mortar.psi_l, T);

				// right element, L2 projection
				L2_projection_to_mortar(temp -> mortar.n_max, temp -> n,
								temp -> index[2], temp -> mortar.l_max, 
								temp -> mortar.a_r, temp -> mortar.b_r,
						 		temp -> solution_int_l, 
								temp -> mortar.psi_r, T);
			
				std::vector<int> index{0, (temp -> mortar.n_max + 1), (temp -> mortar.n_max + 1) * 2};	

				for(int s = 0; s <= temp -> mortar.n_max; ++s){

					// Riemann solver
					Riemann_solver_y(temp -> mortar.psi_l, temp -> mortar.psi_r, 
							temp -> mortar.nflux, -1, index);

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}	

				// L2 project back to element, right element
				L2_projection_to_element(temp -> mortar.n_max, temp -> n, 
							temp -> index[2], temp -> mortar.l_max, 
							temp -> mortar.b_r,
			 				temp -> nflux_l, temp -> mortar.nflux, T);

				// L2 projection from mortar to left element	
				L2_projection_to_element(temp -> mortar.n_max, it_face -> porderx, 
							it_face -> hlevel, temp -> mortar.l_max, 
							temp -> mortar.b_l,
			 				temp -> ghost[n_key], temp -> mortar.nflux, T);

			}
			else if(it_face -> face_type == 'B'){	// phsical boundary

				std::vector<double> solution_ext(size);

				double del_x = ((temp -> xcoords[1]) - (temp -> xcoords[0]));

				std::vector<int> index{0, (porderx + 1), (porderx + 1) * 2};	

				for(int s = 0; s <= porderx; ++s){
					// map to physical plane
					double x = Affine_mapping(nodal::gl_points[porderx][s], (temp -> xcoords[0]), del_x);

					// impose boundary conditions (test) ------------------------------------------------
					//External_state_sin_exact(t, x, (temp -> ycoords[0]), solution_ext, index);
					// ----------------------------------------------------------------------------------

					// impose boundary conditions-------------------------------------------------------
					External_state_Gaussian_exact(t, x, (temp -> ycoords[0]), solution_ext, index);
					// ----------------------------------------------------------------------------------

					// Riemann solver
					Riemann_solver_y(solution_ext, temp -> solution_int_l, 
							temp -> nflux_l, -1, index);
					
					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}
			}
			else{	// mpi boundary

				int n_key = it_face -> key;

				Form_mortar_y(temp, it_face); // allocate space on the mortar
				std::vector<double> T;

				// left element, L2 projection
				L2_projection_to_mortar(temp -> mortar.n_max, it_face -> porderx,
								it_face -> hlevel, temp -> mortar.l_max, 
								temp -> mortar.a_l, temp -> mortar.b_l,
						 		temp -> ghost[n_key], 
								temp -> mortar.psi_l, T);
				// right element, L2 projection
				L2_projection_to_mortar(temp -> mortar.n_max, temp -> n,
								temp -> index[2], temp -> mortar.l_max, 
								temp -> mortar.a_r, temp -> mortar.b_r,
						 		temp -> solution_int_l, 
								temp -> mortar.psi_r, T);

				std::vector<int> index{0, (temp -> mortar.n_max + 1), (temp -> mortar.n_max + 1) * 2};	

				for(int s = 0; s <= temp -> mortar.n_max; ++s){

					// Riemann solver
					Riemann_solver_y(temp -> mortar.psi_l, temp -> mortar.psi_r, 
							temp -> mortar.nflux, -1, index);

					std::transform(index.begin(), index.end(), index.begin(), 
							[](int x){return (x + 1);});		// increment 1
				}

				// L2 project back to element, right element
				L2_projection_to_element(temp -> mortar.n_max, temp -> n, 
							temp -> index[2], temp -> mortar.l_max, 
							temp -> mortar.b_r,
			 				temp -> nflux_l, temp -> mortar.nflux, T);

				// L2 projection from mortar to left element	
				// store remote element's nunerical flux in ghost layer. But first clean up ghost layer.
				std::fill(temp -> ghost[n_key].begin(), temp -> ghost[n_key].end(), 0);

				L2_projection_to_element(temp -> mortar.n_max, it_face -> porderx, 
							it_face -> hlevel, temp -> mortar.l_max, 
							temp -> mortar.b_l,
			 				temp -> ghost[n_key], temp -> mortar.nflux, T);
//if(mpi::rank == 3){
//
//	std::cout.precision(17);
//
//	for(auto& h : temp -> ghost[n_key]){
//
//		std::cout<< std::fixed<< h << "\n";
//	}
//	std::cout<< "\n";
//}


			}
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
	
				// test -------------------------------------------------------------------------------
				//External_state_sin_exact(t, x, (temp -> ycoords[1]), solution_ext, index);
				//------------------------------------------------------------------------------------

				// impose boundary conditions --------------------------------------------------------
				External_state_Gaussian_exact(t, x, (temp -> ycoords[1]), solution_ext, index);
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
