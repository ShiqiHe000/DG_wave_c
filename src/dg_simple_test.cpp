#include "dg_simple_test.h"
#include "dg_unit.h"
#include <vector>
#include "dg_boundary_table.h"
#include "dg_local_storage.h"
#include "dg_boundary_table.h"
#include "dg_cantor_pairing.h"
#include <mpi.h>
#include "dg_elem_length.h"
#include "dg_mpi_table_construct.h"
#include  <iostream>	// test
#include "dg_param.h"	//test


void Simple_test(){

	// init, var = 0 at the beginning 
	

	// tarverse the hash table
	Unit* temp = local::head;

	// construct interface (north to south inferfaces)
	std::vector<int> n_interface(local::local_elem_num);
	std::vector<int> s_interface(local::local_elem_num);
	
	for(int k = 0; k < local::local_elem_num; ++k){

		n_interface[k] = temp -> var;

		temp = temp -> next;
	}


	// exchange info on mpi boundaries

	// start to send info (based on the mpi table)
	


//	if(! north.empty()){
//		for(auto& v : north ){
//			std::cout<< mpi::rank << " " << v.target_rank << "\n";
//		}
//	}



}
