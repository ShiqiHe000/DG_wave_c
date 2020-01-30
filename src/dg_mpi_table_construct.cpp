#include "dg_mpi_table_construct.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_boundary_table.h"
#include "dg_cantor_pairing.h"
#include <algorithm>
#include <vector>

/// @brief 
/// Construct MPI boundary tables.
/// @param north North direction table.
/// @param south South direction table.
// only in x direction
void Construct_mpi_table(std::vector<table_elem>& north, std::vector<table_elem>& south){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		int pre_rank = -1;

		// south
		// interate through face 0
		for(auto& face_s : temp -> facen[0]){


			if(face_s.face_type == 'M' && face_s.rank != pre_rank){	// if mpi boundary and rank changes, record

				south.push_back(table_elem());
				south.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				south.back().target_rank = face_s.rank;
				south.back().coord = temp -> index[1];		// y coord
				south.back().hlevel = temp -> index[2]; 	// hlevel	
			}

		
			pre_rank = face_s.rank;
		
		}
	
		pre_rank = -1;
	
		// north
		// iterate through face 1
		for(auto& face_n : temp -> facen[1]){


			if(face_n.face_type == 'M' && face_n.face_type != pre_rank){	// if mpi boundary, record
				
				north.push_back(table_elem());
				north.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				north.back().target_rank = face_n.rank;
				north.back().coord = temp -> index[1];
				north.back().hlevel = temp -> index[2];	
			}

		
			pre_rank = face_n.rank;
		
		}

		temp = temp -> next;

	}	

		// sort north and south table in the end
		std::sort(south.begin(), south.end(), compare_coord);
		std::sort(north.begin(), north.end(), compare_coord);
}

