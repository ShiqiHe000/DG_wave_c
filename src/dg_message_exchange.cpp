#include <mpi.h>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include <vector>
#include <unordered_map>
#include "dg_boundary_table.h"
#include <algorithm>
#include <cassert>
#include "dg_message_exchange.h" 
#include "dg_param.h"

/// @brief
/// Exchange element interface info for those on the MPI boundaries. In x direction. 
/// North send, south recv.
/// @param
void Exchange_solution(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s,
			std::unordered_map<int, std::vector<mpi_table>>& recver, int face_r, char dir){

	assert(dir == 'x' || dir == 'y' && "Direction can only be 'x' or 'y'.");

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
					// copy solution to a vector to send
					for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){	
					
						std::copy(temp -> solution_int_r[equ].begin(), 
								temp -> solution_int_r[equ].end(), 
								std::back_inserter(send_info));
					}
					int count = send_info.size();
					// tag == sender's key
					MPI_Send(&send_info[0], count, MPI_DOUBLE, target_rank, local_key, MPI_COMM_WORLD);
	
				}
			}
			
		}
		
	}
	
	// recver
	if(dir == 'x'){
		for(auto& v : recver){
			
			int target_rank = v.first;
			auto it_local = v.second.begin();	// point to the members in vactor
	
			for(; it_local != v.second.end(); ++it_local){
	
				int local_key = it_local -> local_key;	// sender's key
				Unit* temp = local::Hash_elem[local_key];
	
				// go to this element, and loop through its facen[face_s]
				for(auto it_face = temp -> facen[face_s].begin();
					it_face != temp -> facen[face_s].end(); ++it_face){
				
					if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ // send
					
						// recv info from target rank
						int recv_size = dg_fun::num_of_equation * (it_face -> pordery + 1);
						int sizei = (it_face -> pordery + 1);
						std::vector<double> recv_info(recv_size);
						
						MPI_Status status;
						MPI_Recv(&recv_info[0], recv_size, MPI_DOUBLE, 
							target_rank, it_face -> key, MPI_COMM_WORLD, &status );
	
						// put it in ghost layer
						auto it1 = recv_info.begin();	// first
						auto it2 = it1 + sizei;	// last
						for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
							temp -> ghost[equ] = std::vector<double>(sizei);
							std::copy(it1, it2, temp -> ghost[equ].begin());
				
							it1 += sizei;
							it2 += sizei;
						}
					}
				}
			}
		}
	}
	else{	// y dir

		for(auto& v : recver){
			
			int target_rank = v.first;
			auto it_local = v.second.begin();	// point to the members in vactor
	
			for(; it_local != v.second.end(); ++it_local){
	
				int local_key = it_local -> local_key;	// sender's key
				Unit* temp = local::Hash_elem[local_key];
	
				// go to this element, and loop through its facen[face_s]
				for(auto it_face = temp -> facen[face_s].begin();
					it_face != temp -> facen[face_s].end(); ++it_face){
				
					if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ // send
					
						// recv info from target rank
						int recv_size = dg_fun::num_of_equation * (it_face -> porderx + 1);
						int sizei = (it_face -> porderx + 1);
						std::vector<double> recv_info(recv_size);
						
						MPI_Status status;
						MPI_Recv(&recv_info[0], recv_size, MPI_DOUBLE, 
							target_rank, it_face -> key, MPI_COMM_WORLD, &status );
	
						// put it in ghost layer
						auto it1 = recv_info.begin();	// first
						auto it2 = it1 + sizei;	// last
						for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
							temp -> ghost[equ] = std::vector<double>(sizei);
							std::copy(it1, it2, temp -> ghost[equ].begin());
				
							it1 += sizei;
							it2 += sizei;
						}
					}
				}
			}
		}
	}

}
