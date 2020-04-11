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
#include "dg_vector_operation.h"
#include "dg_cantor_pairing.h"
#include <iostream>	// test

// forward declaration ----------------------------------------------------------------

void Exchange_solution(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s,
			std::unordered_map<int, std::vector<mpi_table>>& recver, int face_r, char dir);

void Exchange_flux(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s, int face_r);
//--------------------------------------------------------------------------------------



/// @brief
/// Exchange element interface info for those on the MPI boundaries. In x direction. 
/// North send, south recv.
/// @param sender sender's mpi boundary table.
/// @param face_s sender's face number.
/// @param recver recver's mpi boundary table.
/// @param face_r recver's face number.
/// @param dir direction ('x' or 'y').
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
				
					int count = (temp -> solution_int_r).size();
					// tag == sender's key
					MPI_Send(&(temp -> solution_int_r)[0], count, MPI_DOUBLE, 
							target_rank, local_key, MPI_COMM_WORLD);
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
	
				// go to this element, and loop through its facen[face_r]
				for(auto it_face = temp -> facen[face_r].begin();
					it_face != temp -> facen[face_r].end(); ++it_face){
				
					if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ // recv
					
						// recv info from target rank
						int recv_size = dg_fun::num_of_equation * (it_face -> pordery + 1);
						temp -> ghost[it_face -> key] = std::vector<double>(recv_size);
						
						MPI_Status status;
						MPI_Recv(&(temp -> ghost[it_face -> key])[0], recv_size, MPI_DOUBLE, 
							target_rank, it_face -> key, MPI_COMM_WORLD, &status );
//if(mpi::rank == 2){
//
//	std::cout<< "from " << target_rank << " key " << it_face -> key << " size " << recv_size << "\n";
//	
//	for(auto& solu : temp -> ghost[it_face -> key]){
//
//		std::cout<< solu << " ";
//	}
//	std::cout << "\n";
//
//}
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
	
				// go to this element, and loop through its facen[face_r]
				for(auto it_face = temp -> facen[face_r].begin();
					it_face != temp -> facen[face_r].end(); ++it_face){
				
					if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ // send
					
						// recv info from target rank
						int recv_size = dg_fun::num_of_equation * (it_face -> porderx + 1);
						
						temp -> ghost[it_face -> key] = std::vector<double> (recv_size);
						
						MPI_Status status;
						MPI_Recv(&(temp -> ghost[it_face -> key])[0], recv_size, MPI_DOUBLE, 
							target_rank, it_face -> key, MPI_COMM_WORLD, &status );
	
					}
				}
			}
		}
	}

}

/// @brief
/// Get flux on the right interface. 
/// @param sender MPI boundary table of the sender.
/// @param face_s sender's face number.
/// @param face_r receiver's face number. 
void Exchange_flux(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s, int face_r){

	// first send out the fluxes on the south mpi interfaces
	for(auto& v : sender){	// south send

		int target_rank = v.first;
		auto it_local = v.second.begin();	// point to the members in vactor

		for(; it_local != v.second.end(); ++it_local){
	
			int local_key = it_local -> local_key;	// sender's key
			Unit* temp = local::Hash_elem[local_key];

			// go to this element, and loop through its facen[face_s]
			for(auto it_face = temp -> facen[face_s].begin();
				it_face != temp -> facen[face_s].end(); ++it_face){
				
				if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ // send
				
					int n_key = it_face -> key;

					int count = (temp -> ghost[n_key]).size();

					// tag == sender's key
					MPI_Send(&(temp -> ghost[n_key])[0], count, MPI_DOUBLE, 
							target_rank, local_key, MPI_COMM_WORLD);
	
				}
			}
			
		}
		
	}


	// North interface get numerical flux from neighbours
	Unit* temp = local::head;
	for(int k = 0; k < local::local_elem_num; ++k){

		int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		// loop north interface
		for(auto it_face = temp -> facen[face_r].begin(); 
			it_face != temp -> facen[face_r].end(); ++it_face){


			if(it_face -> face_type == 'L'){	// local neighbour

				int n_key = it_face -> key;

				// numerical flux -= neighbours numerical flux
				Vector_minus(local::Hash_elem[n_key] -> ghost[local_key], temp -> nflux_r);

			}
			else if(it_face -> face_type == 'M'){	// remote neighbour
		
				int count = (temp -> nflux_r).size();

				std::vector<double> inter(count); // intermediate vector to store the recv info

				MPI_Status status;

				MPI_Recv(&inter[0], count, MPI_DOUBLE, it_face -> rank,
					 it_face -> key, MPI_COMM_WORLD, &status);	
				
				Vector_minus(inter, temp -> nflux_r);
			}	// 'B' skip
		}


		temp = temp -> next;
	}

}
