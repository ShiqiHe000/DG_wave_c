#include "dg_adapt.h"
#include "dg_hrefinement.h"
#include "dg_local_storage.h"
#include "dg_mpi_table_construct.h"
#include <iostream>	//test
#include <unordered_map>	//test
#include "dg_param.h"	//test

void Adapt(){

	h_refinement();
	
	// x direction
	Construct_mpi_table_x(hrefinement::north, hrefinement::south);
	Update_mpi_boundaries(hrefinement::north, 1, hrefinement::north_accum, hrefinement::south, 0, hrefinement::south_accum);	
	// y direction
	Construct_mpi_table_y(hrefinement::west, hrefinement::east);
	Update_mpi_boundaries(hrefinement::west, 2, hrefinement::west_accum, hrefinement::east, 3, hrefinement::east_accum);	
//if(mpi::rank == 3){
//
//	auto got = local::Hash_elem.find(10);
//
//	if(got != local::Hash_elem.end()){
//
//		std::cout<< local::Hash_elem[10] -> facen[2].front().key << "\n";
//	}
//
//}
}

/// @brief
/// Remove all the elements inside the tables (MPI tables and accumlation table)
void Clear_tables(){

	hrefinement::south_accum.clear();
	hrefinement::north_accum.clear();
	hrefinement::south.clear();
	hrefinement::north.clear();

	hrefinement::west_accum.clear();
	hrefinement::east_accum.clear();
	hrefinement::west.clear();
	hrefinement::east.clear();

}
