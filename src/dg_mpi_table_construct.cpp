#include "dg_mpi_table_construct.h"
#include "dg_local_storage.h"
#include "dg_unit.h"
#include "dg_cantor_pairing.h"
#include <algorithm>
#include <vector>
#include "dg_elem_length.h"
#include <mpi.h>
#include <unordered_map>
#include "dg_param.h"
#include "dg_neighbour_list.h"
#include "dg_boundary_table.h"
#include "dg_put_into_mpi_table.h"
#include "dg_derived_datatype.h"
#include <cassert>	// test
#include <iostream>	// test
#include "dg_write_mpi_table.h"	//test

// forward declaration ------------------------------------------------------------------
void Erase_old_face(std::vector<Unit::Face>::iterator& it_face, std::vector<mpi_table>::iterator& it, 
			int target_rank, int facei);

void Sender_recver(std::unordered_map<int, std::vector<mpi_table>>& south, 
					std::unordered_map<int, std::vector<mpi_table>>& north, int update_dir, 
					std::unordered_map<long long int, std::vector<long long int>>& neighbours_north);

void Update_mpi_boundaries(std::unordered_map<int, std::vector<mpi_table>>& north, int facen,
				std::unordered_map<long long int, std::vector<long long int>>& neighbours_north,
				std::unordered_map<int, std::vector<mpi_table>>& south, int faces, 
				std::unordered_map<long long int, std::vector<long long int>>& neighbours_south);

void Update_hash(std::vector<facen_pack>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& table, 
			int facei, int num, int target_rank, 
			std::unordered_map<long long int, std::vector<long long int>>& neighbours);

void Record_length(int my_hlevel, int n_hlevel, int target_rank, std::unordered_map<int, std::vector<mpi_table>>& my_table);

void Construct_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& north, int face_north, 
				std::unordered_map<long long int, std::vector<long long int>>& neighbours_north, 
				std::unordered_map<int, std::vector<mpi_table>>& south, int face_south, 
				std::unordered_map<long long int, std::vector<long long int>>& neighbours_south);

void Possible_neighbours(Unit* temp, std::unordered_map<long long int, std::vector<long long int>>& neighbours, int facen);

void Clear_tables();

void MPI_table_rebuild(int kt);
//---------------------------------------------------------------------------------------


/// @breif
/// Rebuild MPI table. 
void MPI_table_rebuild(int kt){

	Clear_tables();

	// x direction
	Construct_mpi_table(hrefinement::north, 1, hrefinement::neighbours_north,
				 hrefinement::south, 0, hrefinement::neighbours_south);

	// y direction
	Construct_mpi_table(hrefinement::east, 3, hrefinement::neighbours_east,
				 hrefinement::west, 2, hrefinement::neighbours_west);

}


/// @brief
/// Remove all the elements inside the tables (MPI tables and accumlation table)
void Clear_tables(){

	hrefinement::south.clear();
	hrefinement::north.clear();
	hrefinement::west.clear();
	hrefinement::east.clear();

	hrefinement::neighbours_north.clear();	
	hrefinement::neighbours_south.clear();	
	hrefinement::neighbours_east.clear();	
	hrefinement::neighbours_west.clear();	

}



/// @brief
/// Record the element length that is exposed to a certain neighbour rank. 
/// @param my_level current element's h-refinement level.
/// @param n_hlevel neighbour's h-refinement level.
/// @param target_rank neighbour's rank.
/// @param mpi_table Relevent direction's MPI table. 
void Record_length(int my_hlevel, int n_hlevel, int target_rank, std::unordered_map<int, std::vector<mpi_table>>& my_table){

	if(n_hlevel <= my_hlevel){	// if neighbour is larger or same size
		// length is the full scale of current element
		my_table[target_rank].back().mpi_length = Elem_length(my_hlevel);
	}
	else{	// neighbour is smaller
		my_table[target_rank].back().mpi_length += Elem_length(n_hlevel);

	}

}

/// @brief
/// Form the possible neighbours array based on the face number.
/// @param temp Pointer to the current unit element. 
/// @param neighbours The hash table to store the possible neighbours.
/// @param facen Face number.
void Possible_neighbours(Unit* temp, std::unordered_map<long long int, std::vector<long long int>>& neighbours, int facen){

	long long int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

	if(neighbours.count(local_key) != 0){	// if already record
		
		return;	
	}

	int i = temp -> index[0]; 
	int j = temp -> index[1];
	int k = temp -> index[2];

	// first form the same size neighbour
	if(facen == 0){	// south neighbour

		Neighbours_array_x(i - 1, j, k, local_key, facen, neighbours);

	}
	else if(facen == 1){	// north neighbour
	
		Neighbours_array_x(i + 1, j, k, local_key, facen, neighbours);

	}
	else if(facen == 2){	// west
		Neighbours_array_y(i, j - 1, k, local_key, facen, neighbours);
	}
	else{	// east
	
		Neighbours_array_y(i, j + 1, k, local_key, facen, neighbours);
	}

}


/// @brief 
/// Construct MPI boundary tables. Only in x direction. 
/// @param north MPI boundary table (one direction).
/// @param face_north The face direction of the first MPI boudary table.
/// @param neighbours_north Hash table to store all the possible neighbours (direction north). 
/// @param south MPI boundary table (one direction).
/// @param face_south The face direction of the second MPI boundary table. 
/// @param neighbours_south Hash table to store all the possible neighbours (direction south). 
void Construct_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& north, int face_north, 
				std::unordered_map<long long int, std::vector<long long int>>& neighbours_north, 
				std::unordered_map<int, std::vector<mpi_table>>& south, int face_south, 
				std::unordered_map<long long int, std::vector<long long int>>& neighbours_south){
	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		int pre_rank = -1;
	
		// south
		// interate through face 0
		for(auto it = temp -> facen[face_south].begin(); it != temp -> facen[face_south].end(); ++it){

			if(it -> face_type == 'M'){	// if mpi boundary and rank changes, record

				if(it -> rank != pre_rank){

					Put_in_mpi_table(temp, it, south);
					
					Possible_neighbours(temp, neighbours_south, face_south);
				}

				Record_length(temp -> index[2], it -> hlevel, it -> rank, south);

			}

			pre_rank = it -> rank;
 		
		}

		pre_rank = - 1;

		// north
		// iterate through face 1
		for(auto it = temp -> facen[face_north].begin(); it != temp -> facen[face_north].end(); ++it){

			if(it -> face_type == 'M'){	// if mpi boundary and rank changes, record
				if(it -> rank != pre_rank){

					Put_in_mpi_table(temp, it, north);
					Possible_neighbours(temp, neighbours_north, face_north);
				}

				Record_length(temp -> index[2], it -> hlevel, it -> rank, north);

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
/// @param neighbours_north Hash table to store all the possible neighbours (direction north). 
/// @param south MPI boundary table. 
/// @param faces Face direction of the second table. 
/// @param neighbours_south Hash table to store all the possible neighbours (direction south). 
void Update_mpi_boundaries(std::unordered_map<int, std::vector<mpi_table>>& north, int facen,
				std::unordered_map<long long int, std::vector<long long int>>& neighbours_north,
				std::unordered_map<int, std::vector<mpi_table>>& south, int faces, 
				std::unordered_map<long long int, std::vector<long long int>>& neighbours_south){

	// south send, north recv
	Sender_recver(south, north, facen, neighbours_north);
	// north send, south recv. 
	Sender_recver(north, south, faces, neighbours_south);
}


/// @brief
/// Send and recv info to update MPI boundaries. 
/// @param south sender's MPI boundary table. 
/// @param north recver's MPI boundary table. 
/// @param update_dir update_direaction. Should be the recver's direction.
/// @param neighbours_north Hash table of recver's all possible neighbours. 
void Sender_recver(std::unordered_map<int, std::vector<mpi_table>>& south, 
					std::unordered_map<int, std::vector<mpi_table>>& north, int update_dir, 
					std::unordered_map<long long int, std::vector<long long int>>& neighbours_north){
	int s = south.size();	// number of pairs in the hash table
	int n = north.size();	// number of pairs in the hash table

	if(s > 0){	
		
		for(auto& v : south){	
	
			int target_rank = v.first;
			int num_elem = v.second.size();

			auto it = v.second.begin();

			std::vector<facen_pack> send_info(num_elem);

			// pack the sending info
			for(int k = 0; k < num_elem; ++k){
	
				send_info[k].local_key = it -> local_key;	// key
				send_info[k].hlevel = local::Hash_elem[it -> local_key] -> index[2];	// hlevel
				send_info[k].porderx = local::Hash_elem[it -> local_key] -> n;	// porderx
				send_info[k].pordery = local::Hash_elem[it -> local_key] -> m;	// pordery
				send_info[k].mpi_length = it -> mpi_length;
				send_info[k].Copy_ref(local::Hash_elem[it -> local_key] -> ref_x, 
							local::Hash_elem[it -> local_key] -> ref_y);

				++it;
			}

			MPI_Send(&send_info[0], num_elem, Hash::Facen_type, target_rank, 
					mpi::rank, MPI_COMM_WORLD); 
		}
	}

	//-----------------------------------------------------------------------------------------------------------------------

	// north recv ------------------------------------------------------------------------------------------------------------
	if(n > 0){

		for(auto& v : north){

			MPI_Status status1, status2;		// dummy

			int num{};	// number of elem on the other side

			MPI_Probe(v.first, v.first, MPI_COMM_WORLD, &status1);

			MPI_Get_count(&status1, Hash::Facen_type, &num);
			
			std::vector<facen_pack> recv_info(num);	
	
			MPI_Recv(&recv_info[0], num, Hash::Facen_type, v.first, v.first, MPI_COMM_WORLD, &status2);

			Update_hash(recv_info, north, update_dir, num, v.first, neighbours_north);
		}

	}

	//-----------------------------------------------------------------------------------------------------------------------

}

/// @brief
/// Update facen in Unit hash table.
/// @param recv_info recieved information vector.
/// @param table MPI direction table.
/// @param facei element ith face to be updates
/// @param num recieved element number
/// @param target_rank The rank number of the info sender.
/// @param neighbours possible neighbours hash table. 
void Update_hash(std::vector<facen_pack>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& table, 
			int facei, int num, int target_rank, 
			std::unordered_map<long long int, std::vector<long long int>>& neighbours){
	
	for(auto it = table[target_rank].begin(); it != table[target_rank].end(); ++it){

		// delete old face info			
		auto it_face = local::Hash_elem[it -> local_key] -> facen[facei].begin();
		Erase_old_face(it_face, it, target_rank, facei);

		int l_tol{};

		// find in neighbours list	
		for(auto itn = neighbours[it -> local_key].begin();
			itn != neighbours[it -> local_key].end(); ++itn){

			for(int k = 0; k < num; ++k){
	
				if(recv_info[k].local_key == *(itn)){	// find it!
	
					local::Hash_elem[it -> local_key] -> facen[facei].emplace_back('M', 
								recv_info[k].hlevel, recv_info[k].porderx, 
								recv_info[k].pordery, recv_info[k].local_key,
								target_rank, recv_info[k].ref_x, 
								recv_info[k].ref_y);

					l_tol += recv_info[k].mpi_length;	// record neighbour's length
	
					break;
				}
	
			}
			if(l_tol >= (it -> mpi_length)){break;}
		}

			
	}

}

/// @brief
/// Erase the old mpi boundary info on facei
/// @param it_face iterator on facen.
/// @param it iterator of mpi table.
/// @param facei face direction. 
void Erase_old_face(std::vector<Unit::Face>::iterator& it_face, std::vector<mpi_table>::iterator& it, 
			int target_rank, int facei){

	for(; it_face != local::Hash_elem[it -> local_key] -> facen[facei].end(); ){
		// erase the old info	
		if(it_face -> face_type == 'M' && (it_face -> rank == target_rank)){
			
			it_face = local::Hash_elem[it -> local_key] -> facen[facei].erase(it_face); // it_face move to the next
	
		}
		else{
			++it_face;
		}
	}

}

