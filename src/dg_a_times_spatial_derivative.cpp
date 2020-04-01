#include <vector>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_a_times_spatial_derivative.h"
#include <unordered_map>
#include "dg_spatial_derivative.h"
#include <algorithm>
#include "dg_param.h"
#include "dg_flux_vector.h"
#include "dg_single_index.h"

/// @brief
/// Compute spactial derivetive in x direction. 
/// Here we compute the coefficient matrix A times spatial derivative together. 
void A_times_spatial_derivative_x(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		double del_x = temp -> xcoords[1] - temp -> xcoords[0];

		// allocate space for time derivative
		for(int s = 0; s < dg_fun::num_of_equation; ++s){
			temp -> solution_time_der[s] = std::vector<double> ((temp -> n + 1) * (temp -> m + 1));

		}

		std::vector<int> index_equ{0, temp -> m + 1, (temp -> m + 1) * 2};
		for(int j = 0; j <= (temp -> m); ++j){

			std::unordered_map<int, std::vector<double>> flux_x; // horizontal fluxes <equation_num, xflux>
			std::unordered_map<int, std::vector<double>> flux_der; // <equation_num, flux_der (n + 1)>

			for(int i = 0; i <= (temp -> n); ++i){

				int index = Get_single_index(i, j, (temp -> m + 1));
			
				xflux(temp -> solution, flux_x, index);			
			}
			// flux_der
			Spatial_derivative(temp -> n, flux_x, flux_der, temp, index_equ);

			std::transform(index_equ.begin(), index_equ.end(), index_equ.begin(),
					 [](int x){return (x + 1);});
			
			for(int s = 0; s < dg_fun::num_of_equation; ++s){

				for(int i = 0; i <= (temp -> n); ++i){

					double inter = - (2.0 / del_x) * flux_der[s][i];

					int index = Get_single_index(i, j, (temp -> m + 1));

					(temp -> solution_time_der[s])[index] = inter;

				}
			}

		}
		
		// deallocate solutions on the element boundaries
		(temp -> solution_int_l).clear();
		(temp -> solution_int_r).clear();
		(temp -> nflux_l).clear();
		(temp -> nflux_r).clear();

		temp = temp -> next;
	}
}

/// @brief
/// Compute spactial derivetive in x direction. 
/// Here we compute the coefficient matrix A times spatial derivative together. 
void A_times_spatial_derivative_y(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		double del_y = temp -> ycoords[1] - temp -> ycoords[0];

		std::vector<int> index_equ{0, temp -> n + 1, (temp -> n + 1) * 2};
		for(int i = 0; i <= (temp -> n); ++i){

			std::unordered_map<int, std::vector<double>> flux_y; // horizontal fluxes <equation_num, xflux>
			std::unordered_map<int, std::vector<double>> flux_der; // <equation_num, flux_der (m + 1)>

			for(int j = 0; j <= (temp -> m); ++j){

				int index = Get_single_index(i, j, (temp -> m + 1));
				yflux(temp -> solution, flux_y, index);			
			}
			
			// flux_der
			Spatial_derivative(temp -> m, flux_y, flux_der, temp, index_equ);

			std::transform(index_equ.begin(), index_equ.end(), index_equ.begin(),
					 [](int x){return (x + 1);});
			
			for(int s = 0; s < dg_fun::num_of_equation; ++s){

				for(int j = 0; j <= (temp -> m); ++j){

					double inter = - (2.0 / del_y) * flux_der[s][j];

					int nodei = Get_single_index(i, j, (temp -> m + 1));

	
					(temp -> solution_time_der[s])[nodei] += inter;
				}
			}

		}
		
		// deallocate solutions on the element boundaries
		(temp -> solution_int_l).clear();
		(temp -> solution_int_r).clear();
		(temp -> nflux_l).clear();
		(temp -> nflux_r).clear();

		temp = temp -> next;
	}
}
