#include "dg_adapt.h"
#include "dg_local_storage.h"
#include "dg_mpi_table_construct.h"
#include "dg_amr.h"
#include "dg_error_estimator.h"
#include <iostream>	//test
#include "dg_param.h"	//test
//#include "dg_write_mpi_table.h"


/// @brief
/// Refinement routine. 
/// @param kt kt-th iteration. For testing. 
void Adapt(int kt){

	Flag_elem(kt);

//	Refinement_flag();
	
	hpc_refinement();


//	h_refinement();
//	p_refinement(kt);
	
	// clear the old mpi tables before construct the new ones
	Clear_tables();

	// x direction
	Construct_mpi_table(hrefinement::north, 1, hrefinement::neighbours_north,
				 hrefinement::south, 0, hrefinement::neighbours_south);


	Update_mpi_boundaries(hrefinement::north, 1, hrefinement::neighbours_north,
				hrefinement::south, 0, hrefinement::neighbours_south);	

	// y direction
	Construct_mpi_table(hrefinement::east, 3, hrefinement::neighbours_east,
				 hrefinement::west, 2, hrefinement::neighbours_west);

	Update_mpi_boundaries(hrefinement::east, 3, hrefinement::neighbours_east,
				hrefinement::west, 2, hrefinement::neighbours_west);	


}

