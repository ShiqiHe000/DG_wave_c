#include <mpi.h>
#include "dg_reallocate_elem.h"
#include "dg_local_storage.h"
#include "dg_derived_datatype.h"
#include "dg_param.h"
#include "dg_cantor_pairing.h"
#include <cassert>
#include <iostream> // test

// forward declaration-------------------------------------------------------------------
void Send_pack(std::vector<info_pack>& send_info, std::vector<int>::iterator& it);

void Recv_elem(int source, int tag, std::vector<info_pack>& recv_info, int& count);

void Enlarge_hash(std::vector<info_pack>& recv_info, char dir, int num_recv);

void Face_pack(std::vector<face_pack>& face_info, std::vector<int>& send, int& num);

void Recv_face(int source, int tag, std::vector<face_pack>& recv_face);

void Fill_facen(std::vector<face_pack>& face_info);

void Erase_elem_old(std::vector<int>& send, char dir, int num);
// --------------------------------------------------------------------------------------

/// @brief
/// After built the complete mapping table, now we decide how to reallocate the elements
void Reallocate_elem(){

	int start = LB::elem_accum; 	// first elem global number
	int last = start + local::local_elem_num - 1;	// last elem global number

	int num_pre = LB::Send.pre.size();
	int num_next = LB::Send.next.size();

	MPI_Request request_pre1, request_pre2, request_next1, request_next2;

	if(num_pre > 0){	// something to send
		std::vector<info_pack> send_elem(num_pre);

		auto it = LB::Send.pre.begin();

		std::vector<face_pack> face_info;

		// pack info to send
		Send_pack(send_elem, it);
		int num_n{};
		Face_pack(face_info, LB::Send.pre, num_n);
//	std::cout<<"rank "<< mpi::rank << "\n"	;
		// ready to send 
		MPI_Isend(&send_elem[0], num_pre, Hash::Elem_type, mpi::rank - 1, mpi::rank, MPI_COMM_WORLD, &request_pre1);	// tag = rank
		MPI_Isend(&face_info[0], num_n, Hash::Face_type, mpi::rank - 1, mpi::rank + 1, MPI_COMM_WORLD, &request_pre2);

		Erase_elem_old(LB::Send.pre, 'p', num_pre);
	}
	if(num_next > 0){	

		std::vector<info_pack> send_elem(num_next);

		auto it = LB::Send.next.begin();
		std::vector<face_pack> face_info;

		Send_pack(send_elem, it);
		int num_n{};
		Face_pack(face_info, LB::Send.next, num_n);

		MPI_Isend(&send_elem[0], num_next, Hash::Elem_type, mpi::rank + 1, mpi::rank, MPI_COMM_WORLD, &request_next1);
		MPI_Isend(&face_info[0], num_n, Hash::Face_type, mpi::rank + 1, mpi::rank + 1, MPI_COMM_WORLD, &request_next2);

		Erase_elem_old(LB::Send.next, 'n', num_next);
	}
	
	// recv
	if(mpi::rank == 0){	// first proc
		
		if(LB::proc_mapping_table[1].gnum - 1 > last){	// recv from next

			std::vector<info_pack> recv_info;
			std::vector<face_pack> recv_face;

			int recv_num{};	
			Recv_elem(mpi::rank + 1, mpi::rank + 1, recv_info, recv_num);

			Recv_face(mpi::rank + 1, mpi::rank + 2, recv_face);

			Enlarge_hash(recv_info, 'n', recv_num);
			Fill_facen(recv_face);
		}

	}
	else if(mpi::rank == mpi::num_proc - 1){

		if(LB::proc_mapping_table[mpi::rank].gnum < start){	// recv from pre
			std::vector<info_pack> recv_info;

			std::vector<face_pack> recv_face;

			int recv_num{};	
			Recv_elem(mpi::rank - 1, mpi::rank - 1, recv_info, recv_num);
			
			Recv_face(mpi::rank - 1, mpi::rank, recv_face);

			Enlarge_hash(recv_info, 'p', recv_num);
			Fill_facen(recv_face);

		}

	}
	else{	// proc in between
		if(LB::proc_mapping_table[mpi::rank + 1].gnum - 1 > last){	// recv from next

			std::vector<info_pack> recv_info;
			std::vector<face_pack> recv_face;

			int recv_num{};	
			Recv_elem(mpi::rank + 1, mpi::rank + 1, recv_info, recv_num);
			
			Recv_face(mpi::rank + 1, mpi::rank + 2, recv_face);

			Enlarge_hash(recv_info, 'n', recv_num);
			Fill_facen(recv_face);
		}

		if(LB::proc_mapping_table[mpi::rank].gnum < start){	// recv from pre

			std::vector<info_pack> recv_info;
			std::vector<face_pack> recv_face;
			int recv_num{};	
		
			Recv_elem(mpi::rank - 1, mpi::rank - 1, recv_info, recv_num);

			Recv_face(mpi::rank - 1, mpi::rank, recv_face);

			Enlarge_hash(recv_info, 'p', recv_num);
			Fill_facen(recv_face);
		}
		
	}
	
	// wait
	if(num_pre > 0){
		MPI_Status status;
		MPI_Wait(&request_pre1, &status);
		MPI_Wait(&request_pre2, &status);
	}
	if(num_next > 0){

		MPI_Status status;
		MPI_Wait(&request_next1, &status);
		MPI_Wait(&request_next2, &status);

	}

}

/// @brief
/// After receive element information we create units in Hash table and store them.
/// @parma recv_info The received element info (without facen).
/// @param dir The info coming direction (pre: p, next: 'n').
/// @pram num_recv received element number. 
void Enlarge_hash(std::vector<info_pack>& recv_info, char dir, int num_recv){

	assert((dir == 'p' || dir == 'n') && "The sendign direction can only be 'p' or 'n'.");

	Unit* temp_head = nullptr;	// head pointer to the recved linked list 
	int pre_key;

	for(auto it = recv_info.begin(); it != recv_info.end(); ++it){

		int key = Get_key_fun((*it).index[0], (*it).index[1], (*it).index[2]);

		local::Hash_elem[key] = new Unit();

		if(it == recv_info.begin()){
			temp_head = local::Hash_elem[key];	// temp_head points to the first element. 
		}

		local::Hash_elem[key] -> n = (*it).n;		// for now n == m
		local::Hash_elem[key] -> m = (*it).n;		
		
		for(int i = 0; i < 3; ++i){
			local::Hash_elem[key] -> index[i] = (*it).index[i];
		}

		local::Hash_elem[key] -> status = (*it).status;

		local::Hash_elem[key] -> child_position = (*it).child_position; 

		for(int i = 0; i < 2; ++i){	
			local::Hash_elem[key] -> xcoords[i] = (*it).xcoords[i]; 
			local::Hash_elem[key] -> ycoords[i] = (*it).ycoords[i]; 
		}

		
		local::Hash_elem[key] -> var = (*it).var;
		
		if(it != recv_info.begin()){
			local::Hash_elem[pre_key] -> next = local::Hash_elem[key];

		}

		pre_key = key;
	}

	if(local::local_elem_num == 0){	// if no local_elem

		local::head = temp_head;

	}
	else{

		if(dir == 'p'){		// put the elements at the beginning of the linked list
			
			local::Hash_elem[pre_key] -> next = local::head;
	
			local::head = temp_head;
	
		}	
		else{	// put the element at the end of the linked list
				LB::end -> next = temp_head;
		}

	}

	// update local element number
	local::local_elem_num += num_recv;

}

/// @brief
/// Put the neighbours information inside the hash table. 
/// @param face_info Received face information. 
void Fill_facen(std::vector<face_pack>& face_info){

	for(auto& v : face_info){

		local::Hash_elem[v.owners_key] -> facen[v.facei].push_back(Unit::Face());

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().face_type = v.face_type;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().hlevel = v.hlevel;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().porderx = v.porderx;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().pordery = v.pordery;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().key = v.key;

		local::Hash_elem[v.owners_key] -> facen[v.facei].back().rank = v.rank;

		


	}

}

/// @brief
/// Receive element information from sender. 
/// @param source Sender's rank.
/// @param tag Tag of the message. 
/// @param recv_info Buffer for the message. 
/// @param count number of element received. 
void Recv_elem(int source, int tag, std::vector<info_pack>& recv_info, int& count){

	MPI_Status status1, status2;

	MPI_Probe(source, tag, MPI_COMM_WORLD, &status1);

	MPI_Get_count(&status1, Hash::Elem_type, &count);

	recv_info = std::vector<info_pack>(count);

	MPI_Recv(&recv_info[0], count, Hash::Elem_type, source, tag, MPI_COMM_WORLD, &status2);

}

/// @brief
/// Receive face info.
/// @param source source rank.
/// @param tag message tag.
/// @param recv_face vector to store the message. 
void Recv_face(int source, int tag, std::vector<face_pack>& recv_face){

	MPI_Status status1, status2;

	MPI_Probe(source, tag, MPI_COMM_WORLD, &status1);

	int count;
	MPI_Get_count(&status1, Hash::Face_type, &count);
	
	recv_face = std::vector<face_pack>(count);

	MPI_Recv(&recv_face[0], count, Hash::Face_type, source, tag, MPI_COMM_WORLD, &status2);
}


/// @brief
/// Pack sending information. 
/// @param send_info The sending vector.
/// @param it iterator at the beginning of Sending list (pre or next). 
void Send_pack(std::vector<info_pack>& send_info, std::vector<int>::iterator& it){


	for(auto& v : send_info){
		
		v.n = local::Hash_elem[*it] -> n;

		for(int i = 0; i < 3; ++i){
			v.index[i] = local::Hash_elem[*it] -> index[i];
		}

		v.status = local::Hash_elem[*it] -> status;

		v.child_position = local::Hash_elem[*it] -> child_position;

		for(int i = 0; i < 2; ++i){
		
			v.xcoords[i] = local::Hash_elem[*it] -> xcoords[i];
			v.ycoords[i] = local::Hash_elem[*it] -> ycoords[i];
		}

		v.var = local::Hash_elem[*it] -> var;

		++it;
	}


}

/// @brief
/// Pack face info. 
/// @param face_info Store all the face neighbours in a vector.
/// @param send Sending list. 
/// @param num Number of neighbour. 
void Face_pack(std::vector<face_pack>& face_info, std::vector<int>& send, int& num){

	for(auto& v : send){	// traverse the sending list

		for(int i = 0; i < 4; ++i){	// four faces

			for(auto it = local::Hash_elem[v] -> facen[i].begin(); 
				it != local::Hash_elem[v] -> facen[i].end(); ++it){

				face_info.push_back(face_pack());
				++num;

				face_info.back().owners_key = v; 
				face_info.back().facei = i;
				face_info.back().face_type = it -> face_type;
				face_info.back().hlevel = it -> hlevel;
				face_info.back().porderx = it -> porderx;
				face_info.back().pordery = it -> pordery;
				face_info.back().key = it -> key;
				face_info.back().rank = it -> rank;
			}

		}

	}


}


/// @brief
/// Erase the elements in the sending list from the Hash table
/// @param send Sending list.
/// @param dir Sending direction (pre: 'p', next: 'n').
/// @parma num number of element to be send. 
void Erase_elem_old(std::vector<int>& send, char dir, int num){

	assert((dir == 'p' || dir == 'n') && "The sendign direction can only be 'p' or 'n'.");

	local::local_elem_num -= num;	// remaining elements.

	// take care of linked list
	if(local::local_elem_num == 0){	// nothing left

		local::head = nullptr;
		LB::end = nullptr;

	}
	else{

		if(dir == 'p'){

			int last = send.back();
			local::head = local::Hash_elem[last] -> next;
		}
		else{	// 'n'
			LB::my_rank_last -> next = nullptr;
		}		

	}
	
	// erase elements
	for(auto& v : send){

		local::Hash_elem.erase(v);

	}

}
