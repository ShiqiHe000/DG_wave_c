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
#include "dg_derived_datatype.h"
#include <iostream>	// test

// forward declaration ----------------------------------------------------------------
void Exchange_solution(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s,
			std::unordered_map<int, std::vector<mpi_table>>& recver, int face_r, char dir);

void Exchange_flux(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s, int face_r);

void Pack_send_info(std::vector<double>& send_info, std::vector<double>& source);

void Recv_pairs(int target_rank, int tag, std::vector<neighbour_pair>& recv_pairs);

void Recv_solu_int(int target_rank, int tag, std::vector<double>& recv_solu);

void Exchange_solution_pack(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s,
				std::unordered_map<int, std::vector<mpi_table>>& recver, int face_r, char dir);

void Exchange_flux_pack(std::unordered_map<int, std::vector<mpi_table>>& sender, 
				std::unordered_map<int, std::vector<mpi_table>>& recver, 
				int face_s, int face_r, char dir);
//--------------------------------------------------------------------------------------


/// @brief
/// Non-blocking communication to exchange solutions on the boundaries. 
/// @param sender sender's MPI table. 
/// @param face_s sender's face number. 
/// @param recver recver's MPI table. 
/// @param face_r recver's face number. 
/// @param dir info exchanging direction. Can only by 'x' or 'y'. 
void Exchange_solution_pack(std::unordered_map<int, std::vector<mpi_table>>& sender, int face_s,
				std::unordered_map<int, std::vector<mpi_table>>& recver, int face_r, char dir){

	assert(dir == 'x' || dir == 'y' && "Direction can only be 'x' or 'y'.");

	std::unordered_map<int, std::vector<neighbour_pair>> neighbours;	// neighbours' key, <t_rank, neighbour_pairs>
	std::unordered_map<int, std::vector<double>> solu_send;	// pack all the sending info together

	int num_send = sender.size();	// number of send

	MPI_Request send_request1[num_send];
	MPI_Request send_request2[num_send];

	int isend{};

	for(auto& v : sender){	// North send

		int target_rank = v.first;
		auto it_local = v.second.begin();	// point to the members in vector

		neighbours[target_rank] = std::vector<neighbour_pair>();
		solu_send[target_rank] = std::vector<double>();

		for(; it_local != v.second.end(); ++it_local){
	
			long long int local_key = it_local -> local_key;	// sender's key
			Unit* temp = local::Hash_elem[local_key];

			// go to this element, and loop through its facen[face_s]
			for(auto it_face = temp -> facen[face_s].begin();
				it_face != temp -> facen[face_s].end(); ++it_face){
				
				if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ 
				
					// record neighbours' key
					neighbours[target_rank].emplace_back(local_key, it_face -> key);

					// pack solution
					Pack_send_info(solu_send[target_rank], temp -> solution_int_r);

				}
			}
			
		}

		// send out the info together 
		// adjacent element pairs. tag = self_rank
		int count1 = neighbours[target_rank].size();
		MPI_Isend(&neighbours[target_rank][0], count1, Hash::Adj_pairs, target_rank, mpi::rank, 
				MPI_COMM_WORLD, &send_request1[isend]);

		// solution_int. tag = self_rank + num of proc
		int count2 = solu_send[target_rank].size();
		MPI_Isend(&solu_send[target_rank][0], count2, MPI_DOUBLE, target_rank, 
				mpi::rank + mpi::num_proc, MPI_COMM_WORLD, &send_request2[isend]);

		++isend;
	}

	// recver
	if(dir == 'x'){

		for(auto& v : recver){


			int target_rank = v.first;

			std::vector<neighbour_pair> recv_pairs;	// for adjecent pairs

			// recv neighbour pairs
			Recv_pairs(target_rank, target_rank, recv_pairs);


			// recv solutions on the boundary
			std::vector<double> recv_solu;
			Recv_solu_int(target_rank, target_rank + mpi::num_proc, recv_solu);

			int solui{};

			// match the solution_int with the corresponding elements. 
			for(auto& pair : recv_pairs){

				Unit* temp = local::Hash_elem[pair.recver_key];
	
				// go to this element, and loop through its facen[face_r]
				for(auto it_face = temp -> facen[face_r].begin();
					it_face != temp -> facen[face_r].end(); ++it_face){
				
					if(it_face -> key == pair.sender_key){ // find

						// recv info from target rank
						int recv_size = dg_fun::num_of_equation * (it_face -> pordery + 1);
						temp -> ghost[it_face -> key] = std::vector<double>(recv_size);
					
						for(int s = 0; s < recv_size; ++s){

							temp -> ghost[it_face -> key][s] = recv_solu[solui];

							++solui;
						}
						
						break;	
					}
				}
			}

		}
	}
	else{	// y dir

		for(auto& v : recver){

			int target_rank = v.first;

			std::vector<neighbour_pair> recv_pairs;	// for adjecent pairs

			// recv neighbour pairs
			Recv_pairs(target_rank, target_rank, recv_pairs);

			// recv solutions on the boundary
			std::vector<double> recv_solu;
			Recv_solu_int(target_rank, target_rank + mpi::num_proc, recv_solu);

			int solui{};

			// match the solution_int with the corresponding elements. 
			for(auto& pair : recv_pairs){

				Unit* temp = local::Hash_elem[pair.recver_key];
	
				// go to this element, and loop through its facen[face_r]
				for(auto it_face = temp -> facen[face_r].begin();
					it_face != temp -> facen[face_r].end(); ++it_face){
				
					if(it_face -> key == pair.sender_key){ // find
					
						// recv info from target rank
						int recv_size = dg_fun::num_of_equation * (it_face -> porderx + 1);
						temp -> ghost[it_face -> key] = std::vector<double>(recv_size);
					
						for(int s = 0; s < recv_size; ++s){

							temp -> ghost[it_face -> key][s] = recv_solu[solui];

							++solui;
						}
						
						break;	
					}
				}

			}

		}


	}


	if(num_send > 0){

		MPI_Status status1[num_send];
		MPI_Status status2[num_send];

		MPI_Waitall(num_send, send_request1, status1);
		MPI_Waitall(num_send, send_request2, status2);

	}

}

/// @brief
/// Finction to receive neighbours' solution on the boundary.  
/// @param target_rank Sender's rank number. 
/// @param tag message tag. 
/// @param recv_solu vector to store the message. 
void Recv_solu_int(int target_rank, int tag, std::vector<double>& recv_solu){

	MPI_Status status1, status2;

	MPI_Probe(target_rank, tag, MPI_COMM_WORLD, &status1);

	int count{};

	MPI_Get_count(&status1, MPI_DOUBLE, &count);

	// allocate space
	recv_solu = std::vector<double>(count);

	MPI_Recv(&recv_solu[0], count, MPI_DOUBLE, target_rank, tag, MPI_COMM_WORLD, &status2);	


}


/// @brief
/// Finction to receive neighbour pairs. 
/// @param target_rank Sender's rank number. 
/// @param tag message tag. 
/// @param recv_pairs vector to store the message. 
void Recv_pairs(int target_rank, int tag, std::vector<neighbour_pair>& recv_pairs){

	MPI_Status status1, status2;

	MPI_Probe(target_rank, tag, MPI_COMM_WORLD, &status1);

	int count{};

	MPI_Get_count(&status1, Hash::Adj_pairs, &count);

	// allocate space
	recv_pairs = std::vector<neighbour_pair>(count);

	// recv adjecent pairs
	MPI_Recv(&recv_pairs[0], count, Hash::Adj_pairs, target_rank, tag, MPI_COMM_WORLD, &status2);	


}


/// @brief
/// Pack all the solutions in a large vector and send.
/// @param send_info The vector to store sending info. 
/// @param source The solution to copy from. 
void Pack_send_info(std::vector<double>& send_info, std::vector<double>& source){

	for(auto& v : source){

		send_info.emplace_back(v);

	}
}


/// @brief
/// Exchange element interface info for those on the MPI boundaries. 
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
		auto it_local = v.second.begin();	// point to the members in vector

		for(; it_local != v.second.end(); ++it_local){
	
			long long int local_key = it_local -> local_key;	// sender's key
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
	
				long long int local_key = it_local -> local_key;	// sender's key
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
	
				long long int local_key = it_local -> local_key;	// sender's key
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
/// Get flux on the right interface. Pack fluxes together and send.  
/// @param sender MPI boundary table of the sender.
/// @param face_s sender's face number.
/// @param face_r receiver's face number. 
void Exchange_flux_pack(std::unordered_map<int, std::vector<mpi_table>>& sender, 
				std::unordered_map<int, std::vector<mpi_table>>& recver, 
				int face_s, int face_r, char dir){

	assert(dir == 'x' || dir == 'y' && "Direction can only be 'x' or 'y'.");

	int num_send = sender.size();	// number of send

	MPI_Request send_request1[num_send];
	MPI_Request send_request2[num_send];

	std::unordered_map<int, std::vector<neighbour_pair>> neighbours;	// neighbours' key
	std::unordered_map<int, std::vector<double>> flux_send;	// pack all the sending info together

	int isend{};

	for(auto& v : sender){	// sender

		int target_rank = v.first;
		auto it_local = v.second.begin();	// point to the members in vector

		neighbours[target_rank] = std::vector<neighbour_pair>();

		flux_send[target_rank] = std::vector<double> ();

		for(; it_local != v.second.end(); ++it_local){
	
			long long int local_key = it_local -> local_key;	// sender's key
			Unit* temp = local::Hash_elem[local_key];

			// go to this element, and loop through its facen[face_s]
			for(auto it_face = temp -> facen[face_s].begin();
				it_face != temp -> facen[face_s].end(); ++it_face){
				
				if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ // send
				
					// record neighbours' key
					neighbours[target_rank].emplace_back(local_key, it_face -> key);

					// pack solution
					Pack_send_info(flux_send[target_rank], temp -> ghost[it_face -> key]);

				}
			}
			
		}

		// send out the info together 
		// adjacent element pairs. tag = self_rank
		int count1 = neighbours[target_rank].size();
		MPI_Isend(&neighbours[target_rank][0], count1, Hash::Adj_pairs, target_rank, mpi::rank, 
				MPI_COMM_WORLD, &send_request1[isend]);

		// solution_int. tag = self_rank + num of proc
		int count2 = flux_send[target_rank].size();
		MPI_Isend(&flux_send[target_rank][0], count2, MPI_DOUBLE, target_rank, mpi::rank + mpi::num_proc, 
				MPI_COMM_WORLD, &send_request2[isend]);

		++isend;
	}

	// recv
	if(dir == 'x'){
		for(auto& v : recver){
	
			int target_rank = v.first;
	
			std::vector<neighbour_pair> recv_pairs;	// for adjecent pairs
	
			// recv neighbour pairs
			Recv_pairs(target_rank, target_rank, recv_pairs);
	
			// recv solutions on the boundary
			std::vector<double> recv_flux;
			Recv_solu_int(target_rank, target_rank + mpi::num_proc, recv_flux);
	
			int solui{};
	
			// match the solution_int with the corresponding elements. 
			for(auto& pair : recv_pairs){
	
				Unit* temp = local::Hash_elem[pair.recver_key];
	
				// go to this element, and loop through its facen[face_r]
				for(auto it_face = temp -> facen[face_r].begin();
					it_face != temp -> facen[face_r].end(); ++it_face){
				
					if(it_face -> key == pair.sender_key){ // find
					
						// recv info from target rank
						// recv_size should use local, 
						// neighbour already project the flux to the conforming 
						int recv_size = (temp -> nflux_r).size();
						std::vector<double> inter(recv_size);
		
						for(int s = 0; s < recv_size; ++s){
	
							inter[s] = recv_flux[solui];
	
							++solui;
						}
						
						Vector_minus(inter, temp -> nflux_r);
	
						break;	
					}
				}
	
			}
	
		}
	}
	else{	// y dir

		for(auto& v : recver){
	
			int target_rank = v.first;
	
			std::vector<neighbour_pair> recv_pairs;	// for adjecent pairs
	
			// recv neighbour pairs
			Recv_pairs(target_rank, target_rank, recv_pairs);
	
			// recv solutions on the boundary
			std::vector<double> recv_flux;
			Recv_solu_int(target_rank, target_rank + mpi::num_proc, recv_flux);
	
			int solui{};
	
			// match the solution_int with the corresponding elements. 
			for(auto& pair : recv_pairs){
	
				Unit* temp = local::Hash_elem[pair.recver_key];
	
				// go to this element, and loop through its facen[face_r]
				for(auto it_face = temp -> facen[face_r].begin();
					it_face != temp -> facen[face_r].end(); ++it_face){
				
					if(it_face -> key == pair.sender_key){ // find
					
						// recv info from target rank
						int recv_size = (temp -> nflux_r).size();
						std::vector<double> inter(recv_size);
		
						for(int s = 0; s < recv_size; ++s){
	
							inter[s] = recv_flux[solui];
	
							++solui;
						}
						
						Vector_minus(inter, temp -> nflux_r);
	
						break;	
					}
				}
	
			}
	
		}


	}


	if(num_send > 0){

		MPI_Status status1[num_send];
		MPI_Status status2[num_send];

		MPI_Waitall(num_send, send_request1, status1);
		MPI_Waitall(num_send, send_request2, status2);

	}

	// North interface get numerical flux from neighbours
	Unit* temp = local::head;
	for(int k = 0; k < local::local_elem_num; ++k){

		long long int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		// loop north interface
		for(auto it_face = temp -> facen[face_r].begin(); 
			it_face != temp -> facen[face_r].end(); ++it_face){

			if(it_face -> face_type == 'L'){	// local neighbour

				long long int n_key = it_face -> key;

				// numerical flux -= neighbours numerical flux
				Vector_minus(local::Hash_elem[n_key] -> ghost[local_key], temp -> nflux_r);
			}	// skip 'M' and 'B'
		}


		temp = temp -> next;
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
	
			long long int local_key = it_local -> local_key;	// sender's key
			Unit* temp = local::Hash_elem[local_key];

			// go to this element, and loop through its facen[face_s]
			for(auto it_face = temp -> facen[face_s].begin();
				it_face != temp -> facen[face_s].end(); ++it_face){
				
				if(it_face -> face_type == 'M' && it_face -> rank == target_rank){ // send
				
					long long int n_key = it_face -> key;

					int count = (temp -> ghost[n_key]).size();

					// tag == sender's key
					MPI_Send(&(temp -> ghost[n_key])[0], count, MPI_DOUBLE, 
						target_rank, n_key, MPI_COMM_WORLD);
	
				}
			}
			
		}
		
	}


	// North interface get numerical flux from neighbours
	Unit* temp = local::head;
	for(int k = 0; k < local::local_elem_num; ++k){

		long long int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		// loop north interface
		for(auto it_face = temp -> facen[face_r].begin(); 
			it_face != temp -> facen[face_r].end(); ++it_face){


			if(it_face -> face_type == 'L'){	// local neighbour

				long long int n_key = it_face -> key;

				// numerical flux -= neighbours numerical flux
				Vector_minus(local::Hash_elem[n_key] -> ghost[local_key], temp -> nflux_r);
			}
			else if(it_face -> face_type == 'M'){	// remote neighbour
		
				int count = (temp -> nflux_r).size();

				std::vector<double> inter(count); // intermediate vector to store the recv info

				MPI_Status status;

				MPI_Recv(&inter[0], count, MPI_DOUBLE, it_face -> rank,
					 local_key, MPI_COMM_WORLD, &status);	
				Vector_minus(inter, temp -> nflux_r);
			}	// 'B' skip
		}


		temp = temp -> next;
	}

}
