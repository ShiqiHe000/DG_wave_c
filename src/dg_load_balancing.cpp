#include "dg_load_balancing.h"
#include "dg_reallocate_elem.h"
#include "dg_proc_mapping.h"
#include "dg_local_storage.h"
#include "dg_param.h"	// test
#include <iostream>

// forward declaration--------------------------------
void Clear_mapping_tables();
//----------------------------------------------------

/// @brief
/// Whole procedure of load balancing   
void Load_balancing(int kt){
	Build_mapping_table();

	Update_mpi_boundary(kt);

	Reallocate_elem(kt);
//std::cout<<"rank "<< mpi::rank << " time "<< kt << "\n";
	Clear_mapping_tables();

}

void Clear_mapping_tables(){

	LB::proc_mapping_table.clear();	// proc mapping table

	LB::Send = {};	// Sending list clear up

	LB::elem_accum = 0;

	LB::my_rank_last = nullptr;
	LB::my_rank_first = nullptr;
}
