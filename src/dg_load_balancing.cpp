#include "dg_load_balancing.h"
#include "dg_reallocate_elem.h"
#include "dg_proc_mapping.h"
#include "dg_mpi_table_construct.h"
#include "dg_local_storage.h"
#include "dg_param.h"	// test
#include <iostream>
//#include <TAU.h>	// tau profiling memory

// forward declaration--------------------------------
void Clear_mapping_tables();
void LB_set_back();
void Load_balancing(int kt);
//----------------------------------------------------

/// @brief
/// Whole procedure of load balancing   
void Load_balancing(int kt){
	
//	TAU_PROFILE("load_balancing()", " ", TAU_DEFAULT);
//	TAU_PROFILE_SET_NODE(0);

//	TAU_TRACK_MEMORY_HERE();


	if(LB::first){	// first time evaluate the LB quality
		Build_mapping_table_quality();

	}
	else{
		Build_mapping_table();
	}

	Update_mpi_boundary();

	Reallocate_elem(kt);

	Clear_mapping_tables();

	MPI_table_rebuild();	// rebuild mpi tables

//	TAU_TRACK_MEMORY_HERE();
}

/// @brief
/// Clean up processor mapping tables
void Clear_mapping_tables(){

	LB::proc_mapping_table.clear();	// proc mapping table

	LB::Send = {};	// Sending list clear up

	LB::elem_accum = 0;

	LB::my_rank_last = nullptr;
	LB::my_rank_first = nullptr;

}

/// @brief
/// Set some the parameters in LB algorithm back to original. 
void LB_set_back(){

	LB::first = true;

	LB::high_eff = false;

	LB::load_average = 0.0;

	LB::opt_bottleneck = 0.0;

}
