#include "dg_time_derivative.h"
#include "dg_local_storage.h"
#include <vector>
#include "dg_interface_construct.h"
#include "dg_message_exchange.h"
#include "dg_unit.h"
#include "dg_numerical_flux.h"
#include "dg_a_times_spatial_derivative.h"
#include <iostream>	// test
#include "dg_param.h"	// test

/// @brief
/// Compute the time derivative of all local element.
/// @param t current time step.
void DG_time_der(double t){

	// x direction==================================================================================
	// compute the solutions on the element interfaces
	Unit* temp = local::head;
	for(int k = 0; k < local::local_elem_num; ++k){
	
		Construct_interface_x(temp);

		temp = temp -> next;
	}


	// exchange solution on the mpi boundaries
	Exchange_solution(hrefinement::north, 1, hrefinement::south, 0, 'x');

	// compute the numberical flux
	Numerical_flux_x(t);

	// exchange numerical flux on the mpi boundaries
	Exchange_flux(hrefinement::south, 0, 1);

	// spatial derivative
	A_times_spatial_derivative_x();

//	if(mpi::rank == 0){
//
//		temp = local::head;
////		std::vector<int> index{0, 7, 14};
//	
////		for(int i = 0; i <= 6; ++i){
////
////			std::cout<<index[0]<< " "<< (temp -> solution_int_r[index[0]]) << " "
////				<< index[1] << " "<< (temp -> solution_int_r[index[1]]) << " "
////				<< index[2]<< " "<< (temp -> solution_int_r[index[2]]) << "\n";
////
////			for(auto& v : index){
////
////				++v;
////			}
////
////		}
//
//		std::cout << temp -> solution_time_der[0][0] << "\n";
//
//
//	}
//std::cout<< "rank "<< mpi::rank << "\n";
	//===============================================================================================

	// y direction ==================================================================================
	temp = local::head;
	for(int k = 0; k < local::local_elem_num; ++k){
	
		Construct_interface_y(temp);

		temp = temp -> next;
	}

	// exchange solution on the mpi boundaries
	Exchange_solution(hrefinement::east, 3, hrefinement::west, 2, 'y');

	// compute the numberical flux
	Numerical_flux_y(t);
	
	// exchange numerical flux on the mpi boundaries
	Exchange_flux(hrefinement::west, 2, 3);
	
	// spatial derivative
	A_times_spatial_derivative_y();
	//===============================================================================================
}
