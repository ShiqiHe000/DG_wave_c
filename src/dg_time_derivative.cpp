#include "dg_local_storage.h"
#include <vector>
#include "dg_interface_construct.h"
//#include "dg_unit.h"
//#include "dg_param.h"

/// @brief
/// Compute the time derivative of all local element.
/// @param t current time step.
void DG_time_der(double t){

//	local::solution_int_l = std::vector<std::vector<double>>(local::local_elem_num);
//	local::solution_int_r = std::vector<std::vector<double>>(local::local_elem_num);

	// x direction==================================================================================
	// get the flux on the interfaces

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		Construct_interface_x(temp);

		temp = temp -> next;

	}

	//===============================================================================================
}
