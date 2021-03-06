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
//	Exchange_solution(hrefinement::north, 1, hrefinement::south, 0, 'x');	// exchange element by element
	Exchange_solution_pack(hrefinement::north, 1, hrefinement::south, 0, 'x');	// exchange together

	// compute the numberical flux
	Numerical_flux_x(t);

//if(mpi::rank == 16){
//
//	std::cout<< "compute flux e (x) \n";
//	std::cout<< "exchange flux s (x) \n";
//}
	// exchange numerical flux on the mpi boundaries
//	Exchange_flux(hrefinement::south, 0, 1);
	Exchange_flux_pack(hrefinement::south, hrefinement::north, 0, 1, 'x');

//if(mpi::rank == 16){
//
//	std::cout<< "exchange flux e (x) \n";
//}
	// spatial derivative
	A_times_spatial_derivative_x();

	//===============================================================================================

	// y direction ==================================================================================
	temp = local::head;
	for(int k = 0; k < local::local_elem_num; ++k){
	
		Construct_interface_y(temp);

		temp = temp -> next;
	}

//if(mpi::rank == 0){
//
//	std::cout<< "exchange solu s (y) \n";
//}
	// exchange solution on the mpi boundaries
//	Exchange_solution(hrefinement::east, 3, hrefinement::west, 2, 'y');
	Exchange_solution_pack(hrefinement::east, 3, hrefinement::west, 2, 'y');

	// compute the numberical flux
	Numerical_flux_y(t);
	
	// exchange numerical flux on the mpi boundaries
//	Exchange_flux(hrefinement::west, 2, 3);
	Exchange_flux_pack(hrefinement::west, hrefinement::east, 2, 3, 'y');
	
//if(mpi::rank == 0){
//
//	std::cout<< "exchange solu e (y) \n";
//}
	// spatial derivative
	A_times_spatial_derivative_y();
	//===============================================================================================
}
