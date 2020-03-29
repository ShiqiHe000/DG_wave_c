#include "dg_local_storage.h"
#include <vector>
#include "dg_interface_construct.h"
#include "dg_message_exchange.h"
#include "dg_unit.h"
#include "dg_numerical_flux.h"
//#include "dg_param.h"

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
	
	//===============================================================================================
}
