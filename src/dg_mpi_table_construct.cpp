#include "dg_mpi_table_construct.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_boundary_table.h"
#include "dg_cantor_pairing.h"
#include <algorithm>
#include <vector>
#include "dg_elem_length.h"
#include <mpi.h>
#include <unordered_map>
#include "dg_param.h"
#include <iostream>	// test

// forward declaration ------------------------------------------------------------------
void Erase_old_face(std::vector<Unit::Face>::iterator& it_face, std::vector<mpi_table>::iterator& it, int facei);

void Sender_recver(std::unordered_map<int, std::vector<mpi_table>>& south, 
			std::unordered_map<int, std::vector<mpi_table>>& north, int update_dir);

void Update_mpi_boundaries(std::unordered_map<int, std::vector<mpi_table>>& north, int facen, 
				std::unordered_map<int, std::vector<mpi_table>>& south, int faces);

void Update_hash(std::vector<int>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& table, 
			int facei, int num1, int target_rank);

void Record_length(int my_hlevel, int n_hlevel, int target_rank, std::unordered_map<int, std::vector<mpi_table>>& mpi_table);

void Put_in_mpi_table(Unit* temp, std::vector<Unit::Face>::iterator& facen_it, 
			std::unordered_map<int, std::vector<mpi_table>>& mpi_table)

void Construct_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& north, int face_north
				std::unordered_map<int, std::vector<mpi_table>>& south, int face_south)
//---------------------------------------------------------------------------------------

/// @brief
/// Record the element length that is exposed to a certain neighbour rank. 
/// @param my_level current element's h-refinement level.
/// @param n_hlevel neighbour's h-refinement level.
/// @param target_rank neighbour's rank.
/// @param mpi_table Relevent direction's MPI table. 
void Record_length(int my_hlevel, int n_hlevel, int target_rank, std::unordered_map<int, std::vector<mpi_table>>& mpi_table){

	if(n_hlevel <= my_hlevel){	// if neighbour is larger or same size
		// length is the full scale of current element
		mpi_table[target_rank].back().mpi_length = Elem_length(my_hlevel);
	}
	else{	// neighbour is smaller
		mpi_table[target_rank].back().mpi_length += Elem_length(n_hlevel);

	}

}

/// @brief 
/// Insert the current MPI boundary element to the MPI boundary table. 
/// @param temp Pointer to the current unit. 
/// @param facen_it Iterator on the element face info vector.
/// @param mpi_table The relevent MPI bountary table. 
/// @param target_rank The neighbour's rank. 
void Put_in_mpi_table(Unit* temp, std::vector<Unit::Face>::iterator& facen_it, 
			std::unordered_map<int, std::vector<mpi_table>>& mpi_table){

	if(mpi_table.count(facen_it -> rank) == 0){	// if this rank has not been not record yet

		mpi_table[target_rank] = std::vector<mpi_table>();

		int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		// mpi_length and owners_rank will be recorded later
		mpi_table[target_rank].push_back({local_key, facen_it -> rank, facen_it -> hlevel, 0, 0});
	}
	else{ // this rank is already been record

		int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		mpi_table[target_rank].push_back({local_key, target_rank, facen_it -> hlevel, 0, 0});

	}

}


/// @brief 
/// Construct MPI boundary tables. Only in x direction. 
/// @param north MPI north boundary table.
/// @param south MPI south boundary table.
void Construct_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& north, int face_north
				std::unordered_map<int, std::vector<mpi_table>>& south, int face_south){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		int pre_rank = -1;
	

		// south
		// interate through face 0
		for(auto it = temp -> facen[face_south].begin(); it != temp -> facen[0].end(); ++it){

			if(it -> face_type == 'M'){	// if mpi boundary and rank changes, record

				if(it -> rank != pre_rank){

					Put_in_mpi_table(temp, it, south);
				}

				Record_length(temp -> index[2], it -> hlevel, it -> rank, south);

			}

			pre_rank = it -> rank;
		
		}


		pre_rank = - 1;

		// north
		// iterate through face 1
		for(auto it = temp -> facen[face_north].begin(); it != temp -> facen[1].end(); ++it){

			if(it -> face_type == 'M'){	// if mpi boundary and rank changes, record

				if(it -> rank != pre_rank){

					Put_in_mpi_table(temp, it, north);
				}

				Record_length(temp -> index[2], it -> hlevel, it -> rank, south);

			}

			pre_rank = it -> rank;
		
		}
		
		temp = temp -> next;

	}	

}




/// @brief
/// Updates MPI boundaries
/// @param north MPI boundary table.
/// @param facen Face direction of the first table. 
/// @param south MPI boundary table. 
/// @param faces Face direction of the second table. 
void Update_mpi_boundaries(std::unordered_map<int, std::vector<mpi_table>>& north, int facen, 
				std::unordered_map<int, std::vector<mpi_table>>& south, int faces){

	// south send, north recv
	Sender_recver(south, north, facen);

	// north send, south recv. 
	Sender_recver(north, south, faces);
}


/// @brief
/// Send and recv info to update MPI boundaries. 
/// @param south sender's MPI boundary table. 
/// @param north recver's MPI boundary table. 
/// @param update_dir update_direaction. Should be the recver's direction.
void Sender_recver(std::unordered_map<int, std::vector<mpi_table>>& south, 
					std::unordered_map<int, std::vector<mpi_table>>& north, int update_dir){
	int s = south.size();	// number of pairs in the hash table
	int n = north.size();	// number of pairs in the hash table

	if(s > 0){	
		
		MPI_Request s_request[s];	// for mpi_waitall
		MPI_Status  s_status[s];		// mpi_waitall

		int i{};
		for(auto& v : south){	
	
			int target_rank = v.first;
			int num_elem = v.second.size();
	
			std::vector<int> send_info(num_elem * 5);
			auto it = v.second.begin();
	
			// serialize the struct
			for(int k = 0; k < num_elem; ++k){
	
				send_info[5 * k] = it -> local_key;	// key
				send_info[5 * k + 1] = it -> hlevel;	// hlevel
				send_info[5 * k + 2] = local::Hash_elem[it -> local_key] -> n;	// porderx
				send_info[5 * k + 3] = local::Hash_elem[it -> local_key] -> m;	// pordery
				send_info[5 * k + 4] = it -> mpi_length;
				++it;
			}
	
			MPI_Isend(&send_info[0], num_elem * 5, MPI_INT, target_rank, mpi::rank, MPI_COMM_WORLD, &s_request[i]); 
			++i;
		}
		MPI_Waitall(s, s_request, s_status);
	}

//	if(s > 0){	// there is thing to send
//		MPI_Request s_request[s];	// for mpi_waitall
//		MPI_Status  s_status[s];		// mpi_waitall
//
//		int i{};
//		int j{};
//		for(auto& v : south_accum){	// south send
//
//			std::vector<int> send_info(v.sum * 5);	// key, hlevel, porderx, pordery
//			
//			// serialization the struct
//			for(int k = 0; k < v.sum; ++k){
//				
//				send_info[5 * k] = south[j].local_key;	// key
//				send_info[5 * k + 1] = south[j].hlevel;	// hlevel
//				send_info[5 * k + 2] = local::Hash_elem[south[j].local_key] -> n;	// porderx
//				send_info[5 * k + 3] = local::Hash_elem[south[j].local_key] -> m;	// pordery
//				send_info[5 * k + 4] = south[j].mpi_length;
//				++j;
//			}
//			MPI_Isend(&send_info[0], v.sum * 5, MPI_INT, v.rank, mpi::rank, MPI_COMM_WORLD, &s_request[i]); 
//			++i;
//			
//		}
//		MPI_Waitall(s, s_request, s_status);
//
//	}
	//-----------------------------------------------------------------------------------------------------------------------

	// north recv ------------------------------------------------------------------------------------------------------------
	if(n > 0){

		for(auto& v : north){

			MPI_Status status1, status2;		// dummy

			int num{};	// number of elem on the other side
		
			MPI_Probe(v.first, v.first, MPI_COMM_WORLD, &status1);

			MPI_Get_count(&status1, MPI_INT, &num);
			
			std::vector<int> recv_info(num);	
	
			MPI_Recv(&recv_info[0], num, MPI_INT, v.first, v.first, MPI_COMM_WORLD, &status2);
			
			Update_hash(recv_info, north, update_dir, num, v.first);
		}

	}

//	if(n > 0){
//	
//		auto it = north.begin();	// put the iterator at the begin of the north table
//
//		for(auto& v : north_accum){
//
//			MPI_Status status1, status2;		// dummy
//
//			int num;	// number of elem on the other side
//
//			MPI_Probe(v.rank, v.rank, MPI_COMM_WORLD, &status1);
//
//			MPI_Get_count(&status1, MPI_INT, &num);
//			std::vector<int> recv_info(num);	
//	
//			MPI_Recv(&recv_info[0], num, MPI_INT, v.rank, v.rank, MPI_COMM_WORLD, &status2);
//			Update_hash(recv_info, north, update_dir, num, it, v.rank);	// north recv
//		}
//	}
	//-----------------------------------------------------------------------------------------------------------------------

}

/// @brief
/// Update facen in hash table.
/// @param recv_info recieved information vector.
/// @param table MPI direction table.
/// @param facei element ith face to be updates
/// @param num recieved element number * 4.
/// @param it MPI direction table iterator.
/// @param target_rank The rank number of the info sender.
void Update_hash(std::vector<int>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& table, 
			int facei, int num1, int target_rank){
	
	int l_tot{};

	auto it = table[target_rank].begin();

	int num = num1 / 5;
	for(int k = 0; k < num; ){	// not table but number of recv elem


		int l_local = it -> mpi_length;	// local elem length
		int l_n = recv_info[5 * k + 4];		// neighbour elem length


		if(l_local == l_n){	// if same size

			// delete old face info			
			auto it_face = local::Hash_elem[it -> local_key] -> facen[facei].begin();
			Erase_old_face(it_face, it, facei);

			// type, hlevel, porder, key, rank	
			Unit::Face obj = {'M', recv_info[5 * k + 1], recv_info[5 * k + 2], recv_info[5 * k + 3], 
						recv_info[5 * k], it -> target_rank};	


			local::Hash_elem[it -> local_key] -> facen[facei].emplace(it_face, obj); 

			++k;
			++it;

		}
		else if(l_local < l_n){	// local is smaller
			l_tot = 0;
			// recv_info stall (k), it loop
			// traverse mpi table until face matches
			while(l_tot < l_n && it != table[target_rank].end() ){
				
				// erase old face info
				auto it_face = local::Hash_elem[it -> local_key] -> facen[facei].begin();
				Erase_old_face(it_face, it, facei);

				Unit::Face obj = {'M', recv_info[5 * k + 1], recv_info[5 * k + 2], 
							recv_info[5 * k + 3], recv_info[5 * k], it -> target_rank};

				it_face = local::Hash_elem[it -> local_key] -> facen[facei].emplace(it_face, obj);

				l_tot += it -> mpi_length;
			
				++it;	// go to next local elem
				
				// next element does not face this nieghbour
				if((l_tot + it -> mpi_length) > l_n ){ break; } 
			}

			++k;
			
		}
		else{	// local is larger
			l_tot = 0;
			// erase old face info
			auto it_face = local::Hash_elem[it -> local_key] -> facen[facei].begin();
			Erase_old_face(it_face, it, facei);
			
			// it stall, recv_info loop
			while(k < num && l_tot < l_local){
				Unit::Face obj = {'M', recv_info[5 * k + 1], recv_info[5 * k + 2], 
						recv_info[5 * k + 3], recv_info[5 * k], it -> target_rank};
			
				it_face = local::Hash_elem[it -> local_key] -> facen[facei].emplace(it_face, obj);
				l_tot += recv_info[5 * k + 4];

				++it_face;
				++k;

				if((l_tot + recv_info[5 * k + 4]) > l_local){ break;} 
			}

			++it;
		}
		
	}	

}

/// @brief
/// Erase the old mpi boundary info on facei
/// @param it_face iterator on facen.
/// @param it iterator of mpi table.
/// @param facei face direction. 
void Erase_old_face(std::vector<Unit::Face>::iterator& it_face, std::vector<mpi_table>::iterator& it, int facei){

	int target_rank = it -> target_rank;
	for(; it_face != local::Hash_elem[it -> local_key] -> facen[facei].end(); ){
		// erase the old info	
		if(it_face -> face_type == 'M' && (it_face -> rank == target_rank)){
			
			it_face = local::Hash_elem[it -> local_key] -> facen[facei].erase(it_face); // it_face move to the next
	
			if(it_face -> rank != target_rank){break;}
	
		}
		else{
			++it_face;
		}
	}

}

/// @brief
/// Form rank Accumulate table.
/// @param south MPI south table.
/// @param north MPI north table. 
//void Accum_table(std::unordered_map<int, mpi_table>& south, std::vector<accum_elem>& south_accum){
//
//	for(auto& v : south){
//
//		int rank = v.first;
//		int num = v.second().size();	// number of element
//		south_accum.push_back({rank, num});
//
//	}
//
//
////	if(! south.empty()){	// if not empty
////		
////		south_accum.push_back(accum_elem());
////
////		int rank1 = south.front().target_rank;
////		south_accum.back().rank = rank1;
////
////		for(auto& v : south){
////
////			int rank2 = v.target_rank;
////
////			if(rank2 == rank1){
////
////				south_accum.back().sum += 1;
////				
////			}	
////			else{
////				south_accum.push_back(accum_elem());
////				south_accum.back().rank = rank2;
////				south_accum.back().sum += 1;
////				rank1 = rank2;
////		
////			}		
////			
////		}
////	}
////
//}
//
