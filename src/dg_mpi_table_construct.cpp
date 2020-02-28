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

void Sort_mpi_table(std::vector<table_elem>& north, std::unordered_map<int, std::vector<mpi_table>>& hash_north, 
			std::vector<int>& rank_north);

void Update_hash(std::vector<int>& recv_info, std::unordered_map<int, std::vector<mpi_table>>& table, 
			int facei, int num1, int target_rank);

void Construct_mpi_table_y();
void Construct_mpi_table_x();

//---------------------------------------------------------------------------------------

/// @brief 
/// Construct MPI boundary tables. Only in x direction. 
/// @param north MPI north boundary table.
/// @param south MPI south boundary table.
void Construct_mpi_table_x(){

	Unit* temp = local::head;

	// record all the mpi boundary elements	
	std::vector<table_elem> north;
	std::vector<table_elem> south;

	for(int k = 0; k < local::local_elem_num; ++k){

		int new_added{};
		int pre_rank = -1;
		int local_length = Elem_length(temp -> index[2]);	// current element's length
		double coord_c = (temp -> ycoords[1] + temp -> ycoords[0]) / 2.0;	// central coordinates of the element
		// south
		// interate through face 0
		for(auto& face_s : temp -> facen[0]){


			if(face_s.face_type == 'M' && face_s.rank != pre_rank){	// if mpi boundary and rank changes, record

				south.push_back(table_elem());
				south.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				south.back().target_rank = face_s.rank;
				south.back().coord = coord_c;		// y coord
				south.back().hlevel = temp -> index[2]; 	// hlevel	

				++new_added;
			}
			else if(face_s.face_type == 'L'){

				local_length -= Elem_length(face_s.hlevel);
			
				// central coord of the neighbour
				double coord_n_l = local::Hash_elem[face_s.key] -> ycoords[0];
				double coord_n_r = local::Hash_elem[face_s.key] -> ycoords[1];
				double coord_n_c = ( coord_n_l + coord_n_r) / 2.0;
				if(coord_c > coord_n_c){	// neighbour is on the left side

					coord_c = (temp -> ycoords[1] + coord_n_r) / 2.0;

				}
				else{	// neighbour is on the right side
					coord_c = (temp -> ycoords[0] + coord_n_l) / 2.0;

				}
			}

		
			pre_rank = face_s.rank;
		
		}
		// loop the table reversely to store the mpi_length
		if(new_added > 0){
			for(std::vector<table_elem>::reverse_iterator rit = south.rbegin();
				rit != south.rend(); ++rit){

				rit -> mpi_length = local_length;	// mpi_length
				
				rit -> coord = coord_c;

				--new_added;
	
				if(new_added == 0){break;}	
			}
		}

		pre_rank = -1;
		local_length = Elem_length(temp -> index[2]);	// current element's length
		coord_c = (temp -> ycoords[1] + temp -> ycoords[0]) / 2.0;	// central coordinates of the element

		// north
		// iterate through face 1
		for(auto& face_n : temp -> facen[1]){


			if(face_n.face_type == 'M' && face_n.rank != pre_rank){	// if mpi boundary, record
		
				north.push_back(table_elem());
				north.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				north.back().target_rank = face_n.rank;
				north.back().coord = coord_c;	// x direction
				north.back().hlevel = temp -> index[2];	
				++new_added;

			}
			else if(face_n.face_type == 'L'){

				local_length -= Elem_length(face_n.hlevel);
				
				// central coord of the neighbour
				double coord_n_l = local::Hash_elem[face_n.key] -> ycoords[0];
				double coord_n_r = local::Hash_elem[face_n.key] -> ycoords[1];
				double coord_n_c = ( coord_n_l + coord_n_r) / 2.0;
				if(coord_c > coord_n_c){	// neighbour is on the left side

					coord_c = (temp -> ycoords[1] + coord_n_r) / 2.0;

				}
				else{	// neighbour is on the right side
					coord_c = (temp -> ycoords[0] + coord_n_l) / 2.0;

				}
			}

		
			pre_rank = face_n.rank;
		
		}
		// loop the table reversely to store the mpi_length
		if(new_added > 0){
			for(std::vector<table_elem>::reverse_iterator rit = north.rbegin();
				rit != north.rend(); ++rit){
				rit -> mpi_length = local_length;
				rit -> coord = coord_c;
				--new_added;
	
				if(new_added == 0){break;}	
			}
		}
//if(mpi::rank == 1){
//
//	std::cout << temp -> index[0] << temp -> index[1] << temp -> index[2]<< "\n";
//}
		temp = temp -> next;

	}	

		// sort north and south table in the end
		Sort_mpi_table(south, hrefinement::south, hrefinement::rank_south);
		Sort_mpi_table(north, hrefinement::north, hrefinement::rank_north);
//if(mpi::rank ==3){
//	std::cout << "------------------- \n";
//	std::cout<< "elem_num"<< local::local_elem_num << "\n";
//	for(auto& v : north){
//		
//		std::cout<< "local_key "<< v.local_key << " t_rank"<< v.target_rank << " local_length "<< v.mpi_length<<" ";
//	
//	}
//	std::cout<< "\n";
//	std::cout << "------------------- \n";
//}
}

/// @brief 
/// Construct MPI boundary tables. Only in y direction. 
/// @param north MPI west boundary table.
/// @param south MPI east boundary table.
void Construct_mpi_table_y(){

	Unit* temp = local::head;
	std::vector<table_elem> west;
	std::vector<table_elem> east;

	for(int k = 0; k < local::local_elem_num; ++k){
		int new_added{};
		int pre_rank = -1;
		int local_length = Elem_length(temp -> index[2]);	// current element's length
		double coord_c = (temp -> xcoords[1] + temp -> xcoords[0]) / 2.0;	// central coordinates of the element

		// west
		// interate through face 2
		for(auto& face_s : temp -> facen[2]){


			if(face_s.face_type == 'M' && face_s.rank != pre_rank){	// if mpi boundary and rank changes, record

				west.push_back(table_elem());
				west.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				west.back().target_rank = face_s.rank;
				west.back().coord = coord_c;		// x coord
				west.back().hlevel = temp -> index[2]; 	// hlevel	
				++new_added;
			}
			else if(face_s.face_type == 'L'){

				local_length -= Elem_length(face_s.hlevel);
				
				double coord_n_l = local::Hash_elem[face_s.key] -> xcoords[0];
				double coord_n_r = local::Hash_elem[face_s.key] -> xcoords[1];
				double coord_n_c = ( coord_n_l + coord_n_r) / 2.0;
				if(coord_c > coord_n_c){	// neighbour is on the left side

					coord_c = (temp -> xcoords[1] + coord_n_r) / 2.0;

				}
				else{	// neighbour is on the right side
					coord_c = (temp -> xcoords[0] + coord_n_l) / 2.0;

				}
			}
		
			pre_rank = face_s.rank;
		
		}
		// loop the table reversely to store the mpi_length
		if(new_added > 0){
			for(std::vector<table_elem>::reverse_iterator rit = west.rbegin();
				rit != west.rend(); ++rit){
				
				rit -> mpi_length = local_length;
				rit -> coord = coord_c;
				--new_added;
	
				if(new_added == 0){break;}	
			}
		}

		pre_rank = -1;
		
		local_length = Elem_length(temp -> index[2]);	// current element's length
		coord_c = (temp -> xcoords[0] + temp -> xcoords[1]) / 2.0;	
		// east
		// iterate through face 3
		for(auto& face_n : temp -> facen[3]){


			if(face_n.face_type == 'M' && face_n.rank != pre_rank){	// if mpi boundary, record
		
				east.push_back(table_elem());
				east.back().local_key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
				east.back().target_rank = face_n.rank;
				east.back().coord = coord_c;	// y direction
				east.back().hlevel = temp -> index[2];	
				++new_added;
			}
			else if(face_n.face_type == 'L'){

				local_length -= Elem_length(face_n.hlevel);
				double coord_n_l = local::Hash_elem[face_n.key] -> xcoords[0];
				double coord_n_r = local::Hash_elem[face_n.key] -> xcoords[1];
				double coord_n_c = ( coord_n_l + coord_n_r) / 2.0;
				if(coord_c > coord_n_c){	// neighbour is on the left side

					coord_c = (temp -> xcoords[1] + coord_n_r) / 2.0;

				}
				else{	// neighbour is on the right side
					coord_c = (temp -> xcoords[0] + coord_n_l) / 2.0;

				}
			}

			pre_rank = face_n.rank;
		
		}
		// loop the table reversely to store the mpi_length
		if(new_added > 0){
			for(std::vector<table_elem>::reverse_iterator rit = east.rbegin();
				rit != east.rend(); ++rit){

				rit -> mpi_length = local_length;
	
				rit -> coord = coord_c;

				--new_added;
				if(new_added == 0){break;}	
			}
		}
		temp = temp -> next;

	}	

		// sort north and south table in the end
		Sort_mpi_table(west, hrefinement::west, hrefinement::rank_west);
		Sort_mpi_table(east, hrefinement::east, hrefinement::rank_east);
//if(mpi::rank == 3){
//	std::cout << "--------------------"<< "\n";
//	for(auto& v : west){
//
//		std::cout << "local_key "<<v.local_key<< " rank "<< v.target_rank << " ";
//	}
//	std::cout<< "\n";
//	std::cout << "--------------------"<< "\n";
//}
}


/// @brief sort the MPI_boundary table element in ascending sequence among the same target ranks. 
/// @param north the unsorted MPI boundary table.
/// @param hash_north the MPI boudary table to be formed.
/// @param rank_north the currunt direction's neighbour rank sequnence. 
void Sort_mpi_table(std::vector<table_elem>& north, std::unordered_map<int, std::vector<mpi_table>>& hash_north, 
			std::vector<int>& rank_north){

	if(! north.empty()){
	
		// first sort the element coord in ascending sequence
		std::stable_sort(north.begin(), north.end(), 
				[](const table_elem &left, const table_elem &right){
				return left.coord < right.coord;});

		// then record the coord of first appearance "target_rank" 
		for(auto& v : north){
			if(std::find(rank_north.begin(), rank_north.end(), v.target_rank) == rank_north.end()){	
				// if this rank has not been recored
				rank_north.push_back(v.target_rank);	// record rank
				
				// put it in hash
				hash_north[v.target_rank] = std::vector<mpi_table>();
				// owners_rank will be form in load-balancing
				hash_north[v.target_rank].push_back({v.local_key, v.target_rank, v.hlevel,
									v.mpi_length, 0});
			}
			else{	// this rank already exist

				hash_north[v.target_rank].push_back({v.local_key, v.target_rank, v.hlevel,
									v.mpi_length, 0});

			}

		}

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

//	Accum_table(south, south_accum);
//	Accum_table(north, north_accum);
//if(mpi::rank == 1){
//std::cout << "-------------------------- \n";
//	for(auto& v : south_accum){
//		
//		std::cout << "faces" << faces << "rank " << v.rank << " sum " << v.sum << "\n";
//	}
//std::cout << "-------------------------- \n";
////	for(auto& v : north){
////		std::cout<< "facen" << facen << " l_key " << v.local_key << " t_rank "<< v.target_rank << "\n";
////
////	}
//}
//	int s = south_accum.size();
//	int n = north_accum.size();
	
	// south send, north recv
	Sender_recver(south, north, facen);


	// north send, south recv. 
	Sender_recver(north, south, faces);
//if(mpi::rank == 0){
//
//	std::cout<< " ---------------------------------- \n";
//
//	std::cout<< "check"<< "\n";
//	std::cout<< " ---------------------------------- \n";
//}
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
