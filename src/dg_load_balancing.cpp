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
void Load_balancing(){

	Build_mapping_table();

//if(mpi::rank == 0){
//	std::cout<< "check \n";
//
//}
	Update_mpi_boundary();

	Reallocate_elem();
	Clear_mapping_tables();

}

void Clear_mapping_tables(){

	LB::proc_mapping_table.clear();	// proc mapping table

	LB::Send = {};	// Sending list clear up

	LB::elem_accum = 0;

	Unit* end = nullptr;

	Unit* my_rank_last = nullptr;
}
