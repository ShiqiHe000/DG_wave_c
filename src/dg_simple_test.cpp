#include "dg_simple_test.h"
#include "dg_unit.h"
#include <vector>
#include "dg_boundary_table.h"
#include "dg_local_storage.h"
#include "dg_boundary_table.h"
#include "dg_cantor_pairing.h"
#include <algorithm>	// std::sort
#include <mpi.h>
#include  <iostream>	// test
#include "dg_param.h"	//test

// forward declaration-----------------------------------------------------------------------------
void Construct_mpi_table(std::vector<table_elem>& north, std::vector<table_elem>& south);

void Update_mpi_boundaries(std::vector<table_elem>& north, std::vector<table_elem>& south);

void Accum_table(std::vector<table_elem>& south, std::vector<accum_elem>& south_accum);

void Update_hash(std::vector<int>& recv_info, std::vector<table_elem>& table, 
			int facei, int num, std::vector<table_elem>::iterator& it);

void Sender_recver(int s, int n, std::vector<accum_elem>& south_accum, std::vector<accum_elem>& north_accum, 
			std::vector<table_elem>& south, std::vector<table_elem>& north);
//-------------------------------------------------------------------------------------------------

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
	// first construct table
	std::vector<table_elem> north;
	std::vector<table_elem> south;
	
	Construct_mpi_table(north, south);


	// start to send info (based on the mpi table)
	


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

				south.push_back(table_elem());
				south.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				south.back().target_rank = face_s.rank;
				south.back().coord = temp -> index[1];		// y coord
				south.back().hlevel = temp -> index[2]; 	// hlevel	
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
				north.back().hlevel = temp -> index[2];	
			}

		
		
		}

		temp = temp -> next;

	}	

		// sort north and south table in the end
		std::sort(south.begin(), south.end(), compare_coord);
		std::sort(north.begin(), north.end(), compare_coord);
}


/// @brief
/// Updates MPI boundaries
/// @param north north MPI table
/// @param south south MPI table
void Update_mpi_boundaries(std::vector<table_elem>& north, std::vector<table_elem>& south){

	// accumulates boundaries info
	std::vector<accum_elem> south_accum;
	std::vector<accum_elem> north_accum;

	Accum_table(south, south_accum);
	Accum_table(north, north_accum);

	int s = south_accum.size();
	int n = north_accum.size();
	
	// south send, north recv
	Sender_recver(s, n, south_accum, north_accum, south, north);

	// north send, south recv. 
	Sender_recver(n, s, north_accum, south_accum, north, south);

	// south send -----------------------------------------------------------------------------------------------------------------
//	if(s > 0){	// there is thing to send
//		MPI_Request s_request1[s], s_request2[s];	// for mpi_waitall
//		MPI_Status s_status1[s], s_status2[s];		// mpi_waitall
//
//		// first south send north receive
//		int i{};
//		int j{};
//		for(auto& v : south_accum){	// south send
//			MPI_Isend(&v.sum, 1, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &s_request1[i]);	// tag = local rank
//			
//			std::vector<int> send_info(v.sum * 2);
//			
//			// serialization the struct
//			for(int k = 0; k < v.sum; ++k){
//				
//				send_info[2 * k] = south[j].local_key;	// key
//				send_info[2 * k + 1] = south[j].hlevel;	// hlevel
//				++j;
//	
//			}
//	
//			MPI_Isend(&send_info, v.sum * 2, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &s_request2[i]);
//	
//			++i;
//			
//		}
//	
//		MPI_Waitall(s, s_request1, s_status1);	// ensure all info recved
//		MPI_Waitall(s, s_request2, s_status2);
//
//	}
//	//-----------------------------------------------------------------------------------------------------------------------
//	
//
//	// north recv ------------------------------------------------------------------------------------------------------------
//	if(n > 0){
//		MPI_Status status;		// dummy
//		std::vector<int> recv_info;	// recv: key, hlevel
//		std::vector<table_elem>::iterator it;	// declare an iterator
//		it = north.begin();	// put the iterator at the begin of the north table
//	
//		for(auto& v : north_accum){
//			int num;	// number of elem on the other side
//			MPI_Recv(&num, 1, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &status);
//	
//			recv_info = std::vector<int>(num * 2);
//	
//			MPI_Recv(&recv_info, num * 2, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &status);
//			
//			Update_hash(recv_info, north, 1, num, it);
//		}
//	}
//	//-------------------------------------------------------------------------------------------------------------------------
}

/// @brief
/// Send and recv info to update MPI boundaries. 
/// @param s number of element to send.
/// @param n number of element to recv.
/// @param south_accum sender's accum table.
/// @param north_accum recver's accum table.
/// @param south sender's MPI table.
/// @param north recver's MPI table.
void Sender_recver(int s, int n, std::vector<accum_elem>& south_accum, std::vector<accum_elem>& north_accum, 
			std::vector<table_elem>& south, std::vector<table_elem>& north){


	if(s > 0){	// there is thing to send
		MPI_Request s_request1[s], s_request2[s];	// for mpi_waitall
		MPI_Status s_status1[s], s_status2[s];		// mpi_waitall

		// first south send north receive
		int i{};
		int j{};
		for(auto& v : south_accum){	// south send
			MPI_Isend(&v.sum, 1, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &s_request1[i]);	// tag = local rank
			
			std::vector<int> send_info(v.sum * 2);
			
			// serialization the struct
			for(int k = 0; k < v.sum; ++k){
				
				send_info[2 * k] = south[j].local_key;	// key
				send_info[2 * k + 1] = south[j].hlevel;	// hlevel
				++j;
	
			}
	
			MPI_Isend(&send_info, v.sum * 2, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &s_request2[i]);
	
			++i;
			
		}
	
		MPI_Waitall(s, s_request1, s_status1);	// ensure all info recved
		MPI_Waitall(s, s_request2, s_status2);

	}
	//-----------------------------------------------------------------------------------------------------------------------
	

	// north recv ------------------------------------------------------------------------------------------------------------
	if(n > 0){
		MPI_Status status;		// dummy
		std::vector<int> recv_info;	// recv: key, hlevel
		std::vector<table_elem>::iterator it;	// declare an iterator
		it = north.begin();	// put the iterator at the begin of the north table
	
		for(auto& v : north_accum){
			int num;	// number of elem on the other side
			MPI_Recv(&num, 1, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &status);
	
			recv_info = std::vector<int>(num * 2);
	
			MPI_Recv(&recv_info, num * 2, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &status);
			
			Update_hash(recv_info, north, 1, num, it);
		}
	}
	//-------------------------------------------------------------------------------------------------------------------------

}


/// @brief
/// Update facen in hash table.
/// @param recv_info recieved information vector.
/// @param table MPI direction table.
/// @param facei element ith face to be updates
/// @param num recieved element number.
/// @param it MPI direction table iterator.
void Update_hash(std::vector<int>& recv_info, std::vector<table_elem>& table, 
			int facei, int num, std::vector<table_elem>::iterator& it){
	
			
	// now update the facei neighbour
	auto it_hash = local::Hash_elem[it -> local_key] -> facen[facei].begin();
	int target_rank = it -> target_rank;
	for(; it_hash != local::Hash_elem[it -> local_key] -> facen[facei].end(); ){
		// erase the old info	
		if(it_hash -> face_type == 'M' && (it_hash -> rank == target_rank)){
			
			it_hash = local::Hash_elem[it -> local_key] -> facen[facei].erase(it_hash);

			if(it_hash -> rank != target_rank){break;}

		}
		else{
			++it_hash;
		}
	}
	

	int l_tot{};
	for(int k = 0; k < num; ++k){	// not table but number of recv elem

		int l_local = Elem_length(it -> hlevel);
		l_tot += Elem_length(recv_info[2 * k + 1]);
	
		// type, helvel, porder, key, rank	
		Unit::Face obj = {'M', recv_info[2 * k + 1], 0, recv_info[2 * k], target_rank};	// now we do not have porder info
		
		local::Hash_elem[it -> local_key] -> facen[facei].emplace(it_hash, obj);	// it_hash2 becomes invalid 

		++it_hash;

		if(l_local <= l_tot){	// if neightbour is same size or larger
			++it;
			l_tot = 0;
		}
		else if(k == num - 1){

			++it;
		}

	}	

}



void Accum_table(std::vector<table_elem>& south, std::vector<accum_elem>& south_accum){


	if(! south.empty()){	// if not empty
		
		south_accum.push_back(accum_elem());

		int rank1 = south.front().target_rank;
		south_accum.back().rank = rank1;

		for(auto& v : south){

			int rank2 = v.target_rank;

			if(rank2 == rank1){

				south_accum.back().sum += 1;
				
			}	
			else{
				south_accum.push_back(accum_elem());
				south_accum.back().rank = rank2;
				south_accum.back().sum += 1;
				rank1 = rank2;
		
			}		
			
		}
	}

}
