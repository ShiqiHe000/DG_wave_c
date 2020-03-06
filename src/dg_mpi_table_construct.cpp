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
#include <cmath>
#include <iostream>	// test

// forward declaration ------------------------------------------------------------------
void Erase_old_face(std::vector<Unit::Face>::iterator& it_face, std::vector<mpi_table>::iterator& it, int facei);

void Sender_recver(std::unordered_map<int, std::vector<mpi_table>>& south, 
					std::unordered_map<int, std::vector<mpi_table>>& north, int update_dir, 
					std::unordered_map<int, std::vector<int>>& neighbours_north);

void Update_mpi_boundaries(std::unordered_map<int, std::vector<mpi_table>>& north, int facen,
				std::unordered_map<int, std::vector<int>>& neighbours_north,
				std::unordered_map<int, std::vector<mpi_table>>& south, int faces, 
				std::unordered_map<int, std::vector<int>>& neighbours_south);

void Update_hash(std::vector<int>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& table, 
			int facei, int num1, int target_rank, std::unordered_map<int, std::vector<int>>& neighbours);

void Record_length(int my_hlevel, int n_hlevel, int target_rank, std::unordered_map<int, std::vector<mpi_table>>& mpi_table);

void Put_in_mpi_table(Unit* temp, std::vector<Unit::Face>::iterator& facen_it, 
			std::unordered_map<int, std::vector<mpi_table>>& table);

void Construct_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& north, int face_north, 
				std::unordered_map<int, std::vector<int>>& neighbours_north, 
				std::unordered_map<int, std::vector<mpi_table>>& south, int face_south, 
				std::unordered_map<int, std::vector<int>>& neighbours_south);

void Possible_neighbours(Unit* temp, std::unordered_map<int, std::vector<int>>& neighbours, int facen);

void Total_neighbours(int level_now, int& all, int& my_position);

void Neighbours_array_x(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours);

void Neighbours_array_y(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours);
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
			std::unordered_map<int, std::vector<mpi_table>>& table){

	if(table.count(facen_it -> rank) == 0){	// if this rank has not been not record yet

		table[facen_it -> rank] = std::vector<mpi_table>();

		int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

		// mpi_length and owners_rank will be recorded later
		table[facen_it -> rank].push_back({local_key, facen_it -> rank, temp -> index[2], 0, 0});
	}
	else{ // this rank is already been record

		int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

//if(mpi::rank == 0){
//
//	std::cout<< "rank "<< facen_it -> rank<< "\n";
//	std::cout<< "n_key "<< local_key<< "\n";
//}
		table[facen_it -> rank].push_back({local_key, facen_it -> rank, temp -> index[2], 0, 0});

	}

}


/// @brief
/// Form possible array in h-refinement level decending sequence (x direction).
/// @param i same size neighbour's x integer coordinate.
/// @param j same size neighbour's y integer coordinate.
/// @param k same size neighbour's h-refinement level.
/// @param local_key the key of the current element.
/// @param facen facen face number. 
/// @param neighbours neighbours list. 
void Neighbours_array_x(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours){

	int x{}; 
	
	if(facen == 0){	// south
		x = 1;
	}

	int total_n{};
	int same_size{};

	Total_neighbours(k, total_n, same_size);

	neighbours[local_key] = std::vector<int>(total_n);	// allocate space

	// same size neighbour
	neighbours[local_key][same_size] = Get_key_fun(i, j, k);

	std::vector<int> pre_level{i, j, k};

	int elem = same_size - 1;

	for(int m = k + 1; m <= grid::hlevel_max; ++m){

		int layer_elem = (m - k) * 2;
		
		pre_level[0] = 2 * pre_level[0] + x;
		pre_level[1] = 2 * pre_level[1];
		pre_level[2]++;

		for(int n = layer_elem - 1; n >= 0; --n){
			int key = Get_key_fun(pre_level[0], pre_level[1] + n, pre_level[2]);				
			
			neighbours[local_key][elem] = key;

			elem--;
		}
	}

	pre_level[0] = i / 2;
	pre_level[1] = j / 2;
	pre_level[2] = k - 1;

	elem = same_size + 1;

	for(int m = k - 1; m >= 0; --m){

		int key = Get_key_fun(pre_level[0], pre_level[1], pre_level[2]);

		neighbours[local_key][elem] = key;

		++elem;
	
		pre_level[0] /= 2;
		pre_level[1] /= 2;
		pre_level[2]--;
	}

//if(mpi::rank == 0){
//	std::cout<< "neighbours \n";
//	for(auto& v : neighbours[local_key]){
//
//		std::cout<< v << " "<< "\n";
//	}
//	std::cout<< "\n";
//}
}


/// @brief
/// Form possible array in h-refinement level decending sequence (y direction).
/// @param i same size neighbour's x integer coordinate.
/// @param j same size neighbour's y integer coordinate.
/// @param k same size neighbour's h-refinement level.
/// @param local_key the key of the current element.
/// @param facen facen face number. 
/// @param neighbours neighbours list. 
void Neighbours_array_y(int i, int j, int k, int local_key, int facen, 
			std::unordered_map<int, std::vector<int>>& neighbours){

	int y{};
	if(facen == 2){
		y = 1;
	}
	
	int total_n{};
	int same_size{};

	Total_neighbours(k, total_n, same_size);

	neighbours[local_key] = std::vector<int>(total_n);	// allocate space

	// same size neighbour
	neighbours[local_key][same_size] = Get_key_fun(i, j, k);

	std::vector<int> pre_level{i, j, k};

	int elem = same_size - 1;

	for(int m = k + 1; m <= grid::hlevel_max; ++m){

		int layer_elem = (m - k) * 2;
		
		pre_level[0] = 2 * pre_level[0];
		pre_level[1] = 2 * pre_level[1] + y;
		pre_level[2]++;

		for(int n = layer_elem - 1; n >= 0; --n){
			int key = Get_key_fun(pre_level[0] + n, pre_level[1], pre_level[2]);				
			
			neighbours[local_key][elem] = key;

			elem--;
		}
	}

	pre_level[0] = i / 2;
	pre_level[1] = j / 2;
	pre_level[2] = k - 1;

	elem = same_size + 1;

	for(int m = k - 1; m >= 0; --m){

		int key = Get_key_fun(pre_level[0], pre_level[1], pre_level[2]);

		neighbours[local_key][elem] = key;

		++elem;
	
		pre_level[0] /= 2;
		pre_level[1] /= 2;
		pre_level[2]--;
	}
}


/// @brief
/// Form the possible neighbours array based on the face number.
/// @param temp Pointer to the current unit element. 
/// @param neighbours The hash table to store the possible neighbours.
/// @param facen Face number.
void Possible_neighbours(Unit* temp, std::unordered_map<int, std::vector<int>>& neighbours, int facen){

	int local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

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
	
//if(mpi::rank == 0){
//	std::cout<< "local_key "<< local_key<< "\n";
//}
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
/// Compute the total number of possible neighbours.
/// @param level_now The current element h-refinement level.
/// @param all Total position neighbour number.
/// @param my_position The same size neighbour's number. 
void Total_neighbours(int level_now, int& all, int& my_position){

	for(int i = grid::hlevel_max; i >= 0; --i){

		if(i > level_now){

			all += (int)std::pow(2, i - level_now);

		}
		else if(i == level_now){

			++all;
			my_position = all - 1;
		}
		else{
			++all;
		}
	}

} 

/// @brief 
/// Construct MPI boundary tables. Only in x direction. 
/// @param north MPI boundary table (one direction).
/// @param face_north The face direction of the first MPI boudary table.
/// @param south MPI boundary table (one direction).
/// @param face_south The face direction of the second MPI boundary table. 
/// @param neighbours Hash table to store all the possible neighbours. 
void Construct_mpi_table(std::unordered_map<int, std::vector<mpi_table>>& north, int face_north, 
				std::unordered_map<int, std::vector<int>>& neighbours_north, 
				std::unordered_map<int, std::vector<mpi_table>>& south, int face_south, 
				std::unordered_map<int, std::vector<int>>& neighbours_south){
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
//if(mpi::rank == 0){
//
//	std::cout<< "check \n";
//}
				}

				Record_length(temp -> index[2], it -> hlevel, it -> rank, north);

			}

			pre_rank = it -> rank;
		
		}
		
		temp = temp -> next;

	}	

//if(mpi::rank == 0){
//	std::cout<< "-------------------------------------- \n";
//	for(auto& v: north){
//		std::cout<< "n_rank "<< v.first<< "\n";
//		for(auto& a : v.second){
//
//			std::cout<< a.mpi_length << " ";
//
//		}
//		std::cout<< "\n";
//	}
//	std::cout<< "-------------------------------------- \n";
//}

//if(mpi::rank == 0){
//
//	std::cout<< "-------------------------------------- \n";
//	for(auto& v : neighbours_north){
//
//		std::cout<< "local_key "<< v.first<< "\n";
//		for(auto& a : v.second){
//
//			std::cout<< a << " ";
//
//		}
//
//		std::cout<< "\n";
//
//	}
//	std::cout<< "-------------------------------------- \n";
//}
	
}




/// @brief
/// Updates MPI boundaries
/// @param north MPI boundary table.
/// @param facen Face direction of the first table. 
/// @param south MPI boundary table. 
/// @param faces Face direction of the second table. 
void Update_mpi_boundaries(std::unordered_map<int, std::vector<mpi_table>>& north, int facen,
				std::unordered_map<int, std::vector<int>>& neighbours_north,
				std::unordered_map<int, std::vector<mpi_table>>& south, int faces, 
				std::unordered_map<int, std::vector<int>>& neighbours_south){

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
void Sender_recver(std::unordered_map<int, std::vector<mpi_table>>& south, 
					std::unordered_map<int, std::vector<mpi_table>>& north, int update_dir, 
					std::unordered_map<int, std::vector<int>>& neighbours_north){
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
//if(mpi::rank == 2){
//
//	std::cout<< "local_key "<< send_info[5 * k] << " hlevel "<< send_info[5 * k + 1]<< "\n";
//}
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
			
			Update_hash(recv_info, north, update_dir, num, v.first, neighbours_north);
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
			int facei, int num1, int target_rank, std::unordered_map<int, std::vector<int>>& neighbours){
	
	int num = num1 / 5;
	std::vector<int> recv_bucket(num);
	for(int i = 0; i < num; ++i){
		recv_bucket[i] = i;
	}

	for(auto it = table[target_rank].begin(); it != table[target_rank].end(); ++it){

		// delete old face info			
		auto it_face = local::Hash_elem[it -> local_key] -> facen[facei].begin();
		Erase_old_face(it_face, it, facei);

		int l_tol{};

		// find in neighbours list	
		for(auto itn = neighbours[it -> local_key].begin();
			itn != neighbours[it -> local_key].end(); ++itn){


			for(int k = 0; k < num; ++k){
	
				if(recv_info[5 * k] == *(itn)){	// find it!
	
					// type, hlevel, porder, key, rank	
					Unit::Face obj = {'M', recv_info[5 * k + 1], recv_info[5 * k + 2], 
								recv_info[5 * k + 3], 
								recv_info[5 * k], it -> target_rank};	
			
					it_face = local::Hash_elem[it -> local_key] -> 
							facen[facei].emplace(it_face, obj); 
					++it_face;
				
					l_tol += recv_info[5 * k + 4];	// record neighbour's length
	
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

