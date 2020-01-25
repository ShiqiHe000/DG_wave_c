#include "dg_simple_test.h"
#include "dg_unit.h"
#include <vector>
#include "dg_boundary_table.h"
#include "dg_local_storage.h"
#include "dg_boundary_table.h"
#include "dg_cantor_pairing.h"
#include <algorithm>	// std::sort
#include  <iostream>	// test
#include "dg_param.h"	//test

// forward declaration
void Construct_mpi_table(std::vector<table_elem>& north, std::vector<table_elem>& south);


void Simple_test(){

	// init, var = 0 at the beginning 
	

	// tarverse the hash table
	Unit* temp = local::head;

	// construct interface
	std::vector<int> n_interface(local::local_elem_num);
	std::vector<int> s_interface(local::local_elem_num);
	
	for(int k = 0; k < local::local_elem_num; ++k){

		n_interface[k] = temp -> var;

		temp = temp -> next;
	}


	// exchange info on mpi boundaries
	// first construct table
	std::vector<table_elem> north;
	std::vector<table_elem> south;
	
	Construct_mpi_table(north, south);


	// start to send info

//	if(! north.empty()){
//		for(auto& v : north ){
//			std::cout<< mpi::rank << " " << v.target_rank << "\n";
//		}
//	}
}

// only in x direction
void Construct_mpi_table(std::vector<table_elem>& north, std::vector<table_elem>& south){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){


		// south
		// interate through face 0
		for(auto& face_s : temp -> facen[0]){

			if(face_s.face_type == 'M'){	// if mpi boundary, record

//				int i = temp -> index[0];
//				int j = temp -> index[1];
//				int level = temp -> index[2];

				south.push_back(table_elem());
				south.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				south.back().target_rank = face_s.rank;
				south.back().coord = temp -> index[1];
				
			}

		
		
		}
		
		// north
		// iterate through face 1
		for(auto& face_n : temp -> facen[1]){

			if(face_n.face_type == 'M'){	// if mpi boundary, record
				
				north.push_back(table_elem());
				north.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				north.back().target_rank = face_n.rank;
				north.back().coord = temp -> index[1];
				
			}

		
		
		}

		temp = temp -> next;

	}	

		// sort north and south table in the end
		std::sort(south.begin(), south.end(), compare_coord);
		std::sort(north.begin(), north.end(), compare_coord);
}
