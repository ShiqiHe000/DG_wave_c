#include <mpi.h>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include <vector>
#include <unordered_map>
#include "dg_boundary_table.h"
#include <algorithm>

/// @brief
/// Exchange element interface info for those on the MPI boundaries. In x direction. 
/// North send, south recv.
/// @param
void Exchange_solution_x(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s
			std::unordered_map<int, std::vector<mpi_table>>& recver, int face_r){

	for(auto& v : sender){	// North send

		int target_rank = v.first;
		auto it_local = v.second.begin();	// point to the members in vactor

		for(; it_local != v.second.end(); ++it_local){
	
			int local_key = it_local -> local_key;	// sender's key
			Unit* temp = local::Hash_elem[local_key];

			// go to this element, and loop through its facen[face_s]
			for(auto it_face = temp -> facen[face_s].begin();
				it_face != temp -> facen[face_s].end(); ++it_face){
				
				if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ // send
				
					std::vector<double> send_info;
					for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){	// copy solution to a vector to send
						std::copy(temp -> solution_int_r[equ].begin(), 
								temp -> solution_int_r[equ].end(), 
								std::back_inserter(send_info));
					}
					int count = send_info.size();
					MPI_Send(&send_info[0], count, MPI_DOUBLE, target_rank, int tag, MPI_COMM_WORLD);
	
				}
			}
			
		}
		
	}

}
