#include "dg_proc_mapping.h"
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_param.h"
#include <cmath>
#include <mpi.h>
#include "dg_cantor_pairing.h"
#include <algorithm>
#include "dg_status_table.h"
#include "dg_boundary_table.h"
#include "dg_elem_length.h"
#include "dg_load_struct.h"
#include <unordered_map>
#include <cassert>
#include <iostream>	//test

// forward declaration -----------------------------------------
double Elem_load(int porder);

void Build_mapping_table();

void Update_neighbours();

void Neighbour_change(int facei, int n_key, int my_key, int rank);

void Ownership_one_dir(std::unordered_map<int, std::vector<mpi_table>>& mtable);

void Send_recv_ownership(std::unordered_map<int, std::vector<mpi_table>>& sendo, 
			std::unordered_map<int, std::vector<mpi_table>>& recvo, int facei);

void Update_mpib(std::vector<int>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& otable, 
		int facei, int num1, int target_rank);

void Update_mpi_boundary(int kt);

void Change_face(int num, std::vector<int>& recv_info, std::vector<mpi_table>::iterator& ito, 
			std::vector<Unit::Face>::iterator& it_face);
// ------------------------------------------------------------



/// @brief Calculate the sum of the local computational load.
void Build_mapping_table(){
	
	Unit* temp = local::head;
	
	std::vector<double> lprefix_load(local::local_elem_num);
	int pmapping;

	// local prefix sum of load
	lprefix_load[0] = Elem_load(temp -> n);
	temp = temp -> next;
	for(int k = 1; k < local::local_elem_num; ++k){
		lprefix_load[k] = Elem_load(temp -> n) + lprefix_load[k - 1];
		
		temp = temp -> next;

	}	
	
	double local_load_sum = lprefix_load.back();	// local computational load sum
	double exscan_sum{};	// the load of former processor
	
	// Global prefix sum of load
	MPI_Exscan(&local_load_sum, &exscan_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	// calculate for the average load
	double load_avg{};
	if(mpi::rank == (mpi::num_proc - 1)){ // last proc does the job

		double load_tol = exscan_sum + local_load_sum;

		load_avg = load_tol / mpi::num_proc;

	}
	
	// broadcast average load
	MPI_Bcast(&load_avg, 1, MPI_DOUBLE, mpi::num_proc - 1, MPI_COMM_WORLD);
	
	// Global element number
	MPI_Exscan(&local::local_elem_num, &LB::elem_accum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// form processor mapping table
	temp = local::head;
	LB::my_rank_first = local::head;
	int proc_pre = - 1;
	for(int k = 0; k < local::local_elem_num; ++k){
		
		lprefix_load[k] += exscan_sum;

		pmapping = std::floor((lprefix_load[k] - 0.01) / load_avg);

		assert(pmapping >= 0 && "processor mapping is smaller than 0.");	// check

		// form partial mapping table
		if(pmapping != proc_pre){

			LB::proc_mapping_table.push_back({pmapping, k + LB::elem_accum});
			
			proc_pre = pmapping;
		}	
		
		// form sending_envelope
		if(pmapping < mpi::rank){	// need to be moved to the former proc
			int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

			LB::Send.pre.push_back(key);

			LB::my_rank_first = temp -> next;
		}
		else if(pmapping > mpi::rank){

			int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

			LB::Send.next.push_back(key);
		}
		else{

			LB::my_rank_last = temp;
		}


		temp = temp -> next;
	}

	// send the last element's mapping number to the next rank
	// rank0 ~ rank_max-1 send
	MPI_Request request;
	if(mpi::rank != (mpi::num_proc - 1)){
	
		int last_rank = LB::proc_mapping_table.back().irank;
		MPI_Isend(&last_rank, 1, MPI_INT, mpi::rank + 1, mpi::rank + 1, MPI_COMM_WORLD, &request);// tag == recver's rank

	}
	
	// rank1 ~ rank_max recv
	if(mpi::rank != 0){

		int pre_rank;
		MPI_Status status1;

		MPI_Recv(&pre_rank, 1, MPI_INT, mpi::rank - 1, mpi::rank, MPI_COMM_WORLD, &status1);

		int first_rank = LB::proc_mapping_table.front().irank;

		if(first_rank == pre_rank){	// if equal than erase the first column

			LB::proc_mapping_table.erase(LB::proc_mapping_table.begin());
		}

	}
	
	// wait
	if(mpi::rank != (mpi::num_proc - 1)){
	
		MPI_Status status2;
		MPI_Wait(&request, &status2);
	}

	// call allgather to gather the size of the mapping table on each proc
	int sizet = LB::proc_mapping_table.size();
	std::vector<int> sizea(mpi::num_proc);	// vector to store the sizet
	MPI_Allgather(&sizet, 1, MPI_INT, &sizea[0], 1, MPI_INT, MPI_COMM_WORLD);

	// prepare to build the whole mapping table
	std::vector<int> recvcounts(mpi::num_proc);
	std::vector<int> displs(mpi::num_proc);
	for(int i = 0; i < mpi::num_proc; ++i){
		recvcounts[i] = sizea[i] * 2;
		if(i > 0){
	
			displs[i] = displs[i - 1] + recvcounts[i - 1];
		}

	}
	std::vector<int> senda(sizet * 2);	// send buffer
	for(int i = 0; i < sizet; ++i){
		senda[2 * i] = LB::proc_mapping_table[i].irank;
		senda[2 * i + 1] = LB::proc_mapping_table[i].gnum;

	}

	// form the complete mapping table
	std::vector<int> recv_buff(mpi::num_proc * 2);
	MPI_Allgatherv(&senda[0], sizet * 2, MPI_INT, &recv_buff[0], &recvcounts[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

	// rebuild the mapping table
	LB::proc_mapping_table.clear();	// clear all the elements
	for(int i = 0; i < mpi::num_proc; ++i){
		
		LB::proc_mapping_table.push_back({recv_buff[2 * i], recv_buff[2 * i + 1]});

	}

//if(mpi::rank == 1){
//
//	std::cout<< "--------------------------- \n";
//	for(auto& a : LB::Send.pre){
//		std::cout<< a << " ";
//	}
//	std::cout<<"local elem num "<< local::local_elem_num<< "\n";
//	std::cout<< "\n";
//	std::cout<< "--------------------------- \n";
//}
	// Updates local facen based on teh Sending list
	Update_neighbours();


}

/// @brief
/// Updates the element neighbours based on the Send struct.
void Update_neighbours(){
	
	// first updates the Send pre list --------------------------------------------------------------------------------------------------
	for(auto& key : LB::Send.pre){

		for(int facei = 0; facei < 4; ++facei){

			auto it = local::Hash_elem[key] -> facen[facei].begin();

			// traverse the face
			for(; it != local::Hash_elem[key] -> facen[facei].end(); ++it){

				// skip 'M' and 'B'
				if(it -> face_type == 'L'){

					int n_key = it -> key;	// neighbour's key
					
					// find n_key in pre list
					if(std::find(LB::Send.pre.begin(), LB::Send.pre.end(), n_key) == LB::Send.pre.end()){	// if not find

						// find n_key in the next list
						if(std::find(LB::Send.next.begin(), LB::Send.next.end(), n_key) == LB::Send.next.end()){	// if not find
				
							// This element stays locally, updates
							it -> face_type = 'M';
				//			it -> rank = mpi::rank; //current rank	
							// updates neighbour	
							Neighbour_change(facei, n_key, key, mpi::rank - 1);
							
						}
						else{	// if find
							
							// This element will be sent to next proc
							it -> face_type = 'M';
							it -> rank = mpi::rank + 1;
						
							// updates neighbour
							Neighbour_change(facei, n_key, key, mpi::rank - 1);
						}

					}
					else{	// if find in pre list
						it -> rank = mpi::rank - 1;	// update the rank

					}
						
				}

			}		
		
		}
	} //-------------------------------------------------------------------------------------------------------------


	// updates the next list-----------------------------------------------------------------------------------------
	for(auto& key : LB::Send.next){

		for(int facei = 0; facei < 4; ++facei){
			
			auto it = local::Hash_elem[key] -> facen[facei].begin();

			// traverse the face
			for(; it != local::Hash_elem[key] -> facen[facei].end(); ++it){

				// skip 'M' and 'B'
				if(it -> face_type == 'L'){

					int n_key = it -> key;	// neighbour's key

					// find n_key in the next list
					if(std::find(LB::Send.next.begin(), LB::Send.next.end(), n_key) == LB::Send.next.end()){	// if not find
			
						// This element stays locally, updates
						it -> face_type = 'M';
						it -> rank = mpi::rank; //current rank	
						// updates neighbour	
						Neighbour_change(facei, n_key, key, mpi::rank + 1);
						
					}
					else{
						it -> rank = mpi::rank + 1;
					}

				}
			}		

		}

	}
	// --------------------------------------------------------------------------------------------------------------

}

/// @brief
/// Change neighbour's corresponding face. The neighbour element should store locally. 
/// @param facei current face direction. 
/// @param n_key neighbour's key. 
/// @param my_key current element's key. 
/// @param rank Rank number that the current element will reside on. 
void Neighbour_change(int facei, int n_key, int my_key, int rank){

	// neighbour's face direction.
	int oface = Opposite_dir(facei);

	for(auto it = local::Hash_elem[n_key] -> facen[oface].begin(); it != local::Hash_elem[n_key] -> facen[oface].end(); ++it){

		if(it -> key == my_key){

			it -> face_type = 'M';
			it -> rank = rank;

			break;
		}

	}

}


/// @brief
/// Fill in ownership in MPI table in one direction.
/// @param mtable MPI boundary table of the corresponding direction. 
void Ownership_one_dir(std::unordered_map<int, std::vector<mpi_table>>& mtable){


	// keys are inherited from MPI boundary tables
	for(auto& v1 : mtable){

		for(auto& v2 : v1.second){
			if(std::find(LB::Send.pre.begin(), LB::Send.pre.end(), v2.local_key) != LB::Send.pre.end()){ // if find in pre list
				v2.owners_rank = mpi::rank - 1;
			}
			else if(std::find(LB::Send.next.begin(), LB::Send.next.end(), v2.local_key) != LB::Send.next.end()){	// if find in next list
	
				v2.owners_rank = mpi::rank + 1;
	
			}
			else{ // not inside the sending list, record directly
				v2.owners_rank = mpi::rank;
			}
		}

	}

}

/// @brief
/// Updates the MPI boundaries before repartitioning. 
void Update_mpi_boundary(int kt){

	// form the element future ownership----------------------------------------------
	Ownership_one_dir(hrefinement::north);
	Ownership_one_dir(hrefinement::south);
	
	Ownership_one_dir(hrefinement::west);
	Ownership_one_dir(hrefinement::east);
	//-------------------------------------------------------------------------

	// x direction-----------------------------------------------------------------------------------
	
	// north send and south recv
	Send_recv_ownership(hrefinement::north, hrefinement::south, 0);
	// south send and north recv
	Send_recv_ownership(hrefinement::south, hrefinement::north, 1);
	//-----------------------------------------------------------------------------------------------

	// y direction-----------------------------------------------------------------------------------
	
	// west send and east recv
	Send_recv_ownership(hrefinement::west, hrefinement::east, 3);

	// east send and west recv
	Send_recv_ownership(hrefinement::east, hrefinement::west, 2);
	//-----------------------------------------------------------------------------------------------
//std::cout<< "rank "<< mpi::rank<< " kt "<< kt<< "\n";
}

/// @brief
/// Send and recv MPI boundary table to updates the MPI boundaries. 
/// @param sendo Sender's MPI boundary table.
/// @param recvo Recver's MPI boundary table.
/// @param send_accum Sender's accumulation table. 
/// @param recv_accum Receiver's accumulation table. 
/// @param facei the face direction to be updated. 
void Send_recv_ownership(std::unordered_map<int, std::vector<mpi_table>>& sendo, 
			std::unordered_map<int, std::vector<mpi_table>>& recvo, int facei){
	
	int size_s = sendo.size();
	int size_r = recvo.size();
	
	// sender
	if(size_s > 0){	// there is something to send
		MPI_Request s_request[size_s];
		MPI_Status s_status[size_s];

		int i{};
		for(auto& v : sendo){

			int num_elem = v.second.size();	
			std::vector<int> send_info(num_elem * 2);
			auto it = v.second.begin();
	
			// serialization the struct
			for(int k = 0; k < num_elem; ++k){

				send_info[2 * k] = it -> local_key;	
				send_info[2 * k + 1] = it -> owners_rank;
				++it;
			}
assert(v.first >= 0 && v.first < 4 && "target_rank wrong" );
			MPI_Isend(&send_info[0], num_elem * 2, MPI_INT, v.first, mpi::rank, MPI_COMM_WORLD, &s_request[i]);

			++i;

		}

		MPI_Waitall(size_s, s_request, s_status);
	}


	// recver
	if(size_r > 0){

		for(auto& v : recvo){
			
			MPI_Status status1, status2;

			int num{};
assert(v.first >= 0 && v.first <= 3 && "rank wrong send_recv_ownership");
			MPI_Probe(v.first, v.first, MPI_COMM_WORLD, &status1);

			MPI_Get_count(&status1, MPI_INT, &num);

			std::vector<int> recv_info(num);

			MPI_Recv(&recv_info[0], num, MPI_INT, v.first, v.first, MPI_COMM_WORLD, &status2);
//if(mpi::rank == 1){
//	for(int m = 0; m < num / 4; ++m){
//		std::cout<< "key "<< recv_info[4 * m]<< " owner "<< recv_info[4 * m + 3]<< " ";
//	}
//	std::cout<< "\n";
//}			
			Update_mpib(recv_info, recvo, facei, num, v.first);
//std::cout<< "rank "<< mpi::rank << "\n";
		}

	}

}

/// @brief
/// Update MPI boundaries based on the ownership table.
/// @param recv_info Received infomation.
/// @param otable MPI boundary table.
/// @param ito iterator of the MPI boundary table.
/// @param num1 number of element received * 2.
void Update_mpib(std::vector<int>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& otable, 
		int facei, int num1, int target_rank){

	int num = num1 / 2;

	for(auto ito = otable[target_rank].begin(); ito != otable[target_rank].end(); ++ito){

		for(auto it_face = local::Hash_elem[ito -> local_key] -> facen[facei].begin();
			it_face != local::Hash_elem[ito -> local_key] -> facen[facei].end(); ++it_face){

			if(it_face -> face_type == 'M' && it_face -> rank == target_rank){

				Change_face(num, recv_info, ito, it_face);
			}

		}

	}

}

/// @brief
/// Compare the ownership of neighbour element and the stored rank to decide whether to update the facen.
/// @param k k-th element inside the received info.
/// @param recv_info Received information.
/// @param ito interator of the updating MPI boundary table.
/// @param it_face iterator of facen at the corresponding position. 
void Change_face(int num, std::vector<int>& recv_info, std::vector<mpi_table>::iterator& ito, 
			std::vector<Unit::Face>::iterator& it_face){

	for(int k = 0; k < num; ++k){	// loop to find the neighbour

		if(recv_info[2 * k] == it_face -> key){	// find it

			if(recv_info[2 * k + 1] == (ito -> owners_rank)){	// if will be in the same rank
				
				// change 'M' to 'L'
				it_face -> face_type = 'L';
				it_face -> rank = ito -> owners_rank;
			}
			else{	// not in the same rank
				int rank_old = it_face -> rank;
		
				if(rank_old != recv_info[2 * k + 1]){	// if this element will to assign to a new rank
					// change to target rank
					it_face -> rank = recv_info[2 * k + 1];
		
				}	// else no change
			}

		}
	}
}



/// @brief 
/// Element computational load. The load on each element due to fluid computations is O(N**4),
/// where N is the number of grid points along one direction. The load is normalized between 1 and 2. 
/// @param porder input the polynomial order of current element. (Now we assume polynomial
/// order are identical on two direction.)
double Elem_load(int porder){

	static const double load_min = std::pow((double)(grid::nmin + 1), 4);
	static const double load_max = std::pow((double)(grid::nmax + 1), 4);
	
	double load = (std::pow(porder + 1, 4) - load_min) / (load_max - load_min) + 1.0;

	return load;

}
