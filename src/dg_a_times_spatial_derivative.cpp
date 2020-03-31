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
#include <iostream>

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
//if(mpi::rank == 0){
////	std::cout<< "index "<< index<< "\n";
////	std::cout<< "i "<< i << " j "<< j << " equ "<< s << " "<< temp -> solution_time_der[s][index] << "\n";
//
//	std::cout<< "inside "<< temp -> solution_time_der[0][0] << "\n";
//	std::cout<< " ------------------------------------- \n";
//}

				}
			}

		}
		
//if(mpi::rank == 0){
////	std::cout<< "index "<< index<< "\n";
//	std::cout<< temp -> solution_time_der[0][0] << "\n";
//
//}
		// deallocate solutions on the element boundaries
		(temp -> solution_int_l).clear();
		(temp -> solution_int_r).clear();
		(temp -> nflux_l).clear();
		(temp -> nflux_r).clear();

		temp = temp -> next;
	}
//if(mpi::rank == 0){
//	std::cout <<"-------------------------------------------- \n";
//	temp = local::head;
//	std::cout<< "outside "<<temp -> solution_time_der[0][0] << "\n";
//
////	std::cout<< "i "<< i << " j "<< j << " equ "<< s << " "<< temp -> solution_time_der[s][index] << "\n";
//
//}
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
			
//if(mpi::rank == 0){
//	std::cout << " y below------------------------------------------------ \n";
//
//}
			// flux_der
			Spatial_derivative(temp -> m, flux_y, flux_der, temp, index_equ);

//if(mpi::rank == 0){
//
////	if(j == 0){
//		for(int h = 0; h <= 6; ++h){
//			std::cout<<i <<" "<< h << " "<< flux_y[0][h] << " "<< flux_y[1][h]<<" "<<flux_y[2][h] << "\n";
//		}
////	} 
//
//}
			std::transform(index_equ.begin(), index_equ.end(), index_equ.begin(),
					 [](int x){return (x + 1);});
			
			for(int s = 0; s < dg_fun::num_of_equation; ++s){

				for(int j = 0; j <= (temp -> m); ++j){

					double inter = - (2.0 / del_y) * flux_der[s][j];

					int nodei = Get_single_index(i, j, (temp -> m + 1));

//if(mpi::rank == 0){
//
//	std::cout<< "i "<< i << " j "<< j << " s "<< s << " "<< temp -> solution_time_der[s][nodei]<< " inter "<< inter<< "\n";
//}
	
					(temp -> solution_time_der[s])[nodei] += inter;
//if(mpi::rank == 0){
//
////	for(int h = 0; h <= 6; ++h){
//		std::cout<<i <<" "<< j << " "<< s << " "<< temp -> solution_time_der[s][nodei]<< " flux_der "
//			<< flux_der[s][j]<< "\n";
////	}
//
//}
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
