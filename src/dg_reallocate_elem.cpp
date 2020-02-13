#include <mpi.h>
#include "dg_reallocate_elem.h"
#include "dg_local_storage.h"
#include "dg_derived_datatype.h"
#include "dg_param.h"
#include "dg_cantor_pairing.h"

// forward declaration-------------------------------------------------------------------
void Send_pack(std::vector<info_pack>& send_info, std::vector<int>::iterator& it);

void Recv(int source, int tag, std::vector<info_pack>& recv_info);

void Enlarge_hash(std::vector<info_pack>& recv_info, char dir);
// --------------------------------------------------------------------------------------

/// @brief
/// After built the complete mapping table, now we decide how to reallocate the elements
void Reallocate_elem(){

	int start = LB::elem_accum; 	// first elem global number
	int last = start + local::local_elem_num - 1;	// last elem global number
//std::cout<< "rank " << mpi::rank << " start "<< start << " last " << last << "\n";

	int num_pre = LB::Send.pre.size();
	int num_next = LB::Send.next.size();
//std::cout<< "rank " << mpi::rank << "pre " << num_pre << " next " << num_next << "\n";

	MPI_Request request_pre, request_next;

	if(num_pre > 0){	// something to send
		std::vector<info_pack> send_elem(num_pre);

		auto it = LB::Send.pre.begin();

		// pack info to send
		Send_pack(send_elem, it);
	
		// ready to send 
		MPI_Isend(&send_elem[0], num_pre, Hash::Elem_type, mpi::rank - 1, mpi::rank, MPI_COMM_WORLD, &request_pre);	// tag = rank
			
	}
	if(num_next > 0){	

		std::vector<info_pack> send_elem(num_next);

		auto it = LB::Send.next.begin();

		Send_pack(send_elem, it);

		MPI_Isend(&send_elem[0], num_next, Hash::Elem_type, mpi::rank + 1, mpi::rank, MPI_COMM_WORLD, &request_next);

	}

	// recv
	if(mpi::rank == 0){	// first proc
		
		if(LB::proc_mapping_table[1].gnum - 1 > last){	// recv from next

			std::vector<info_pack> recv_info;

			Recv(mpi::rank - 1, mpi::rank - 1, recv_info);

			Enlarge_hash(recv_info, 'n');
		}

	}
	else if(mpi::rank == mpi::num_proc - 1){

		if(LB::proc_mapping_table[mpi::rank].gnum < start){	// recv from pre

			std::vector<info_pack> recv_info;

			Recv(mpi::rank - 1, mpi::rank - 1, recv_info);
			
			Enlarge_hash(recv_info, 'p');
//if(mpi::rank == 3){
//
//	for(auto& n : recv_info){
//		std::cout << "n: " << n.n << "\n";
//		std::cout << "index: " << n.index[0] << n.index[1] << n.index[2] << "\n";
//
//		std::cout<< "status " << n.status << "\n";
//
//		std::cout<< "child " << n.child_position << "\n";
//
//		std::cout << "coordx " << n.xcoords[0] << " " << n.xcoords[1] << "\n";
//		std::cout << "coordy " << n.ycoords[0] << " " << n.ycoords[1] << "\n";
//
//		std::cout << "var " << n.var << "\n";
//		std::cout << "---------------------------------------------------" << "\n";
//
//	}
//
//}

		}

	}
	else{	// proc in between
		if(LB::proc_mapping_table[mpi::rank + 1].gnum - 1 > last){	// recv from next

			std::vector<info_pack> recv_info;

			Recv(mpi::rank + 1, mpi::rank + 1, recv_info);
			
			Enlarge_hash(recv_info, 'n');
		}

		if(LB::proc_mapping_table[mpi::rank].gnum < start){	// recv from pre

			std::vector<info_pack> recv_info;

			Recv(mpi::rank - 1, mpi::rank - 1, recv_info);

			Enlarge_hash(recv_info, 'p');
		}
		
	}
	
	// wait
	if(num_pre > 0){
		MPI_Status status;
		MPI_Wait(&request_pre, &status);
	}

	if(num_next > 0){

		MPI_Status status;
		MPI_Wait(&request_next, &status);

	}

}

/// @brief
/// After receive element information we create units in Hash table and store them.
/// @parma recv_info The received element info (without facen).
/// @param dir The info coming direction (pre: p, next: 'n').
void Enlarge_hash(std::vector<info_pack>& recv_info, char dir){

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


	if(dir == 'p'){		// put the elements at the beginning of the linked list
		
		local::Hash_elem[pre_key] -> next = local::head;

		local::head = temp_head;

	}	
	else if(dir == 'n'){	// put the element at the end of the linked list
		if(local::local_elem_num == 0){
			local::head = temp_head;
		}
		else{	
			LB::end -> next = temp_head;
		}
	}

}

/// @brief
/// Receive element information from sender. 
/// @param source Sender's rank.
/// @param tag Tag of the message. 
/// @param recv_info Buffer for the message. 
void Recv(int source, int tag, std::vector<info_pack>& recv_info){

	MPI_Status status1, status2;
	int count;

	MPI_Probe(source, tag, MPI_COMM_WORLD, &status1);

	MPI_Get_count(&status1, Hash::Elem_type, &count);

	recv_info = std::vector<info_pack>(count);

	MPI_Recv(&recv_info[0], count, Hash::Elem_type, source, tag, MPI_COMM_WORLD, &status2);

}

/// @brief
/// Pack sending information. 
/// @param send_info The sending vector.
/// @param it iterator at the beginning of Sending list (pre or next). 
/// @param dir sending direction (pre: 'p', next: 'n')
void Send_pack(std::vector<info_pack>& send_info, std::vector<int>::iterator& it, char dir){

	if(dir == 'p'){


	}


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

		// erase the element after we recorded it
		local::Hash.erase(*it);
		
		--local::local_elem_num; 

		++it;
	}

	if(local::local_elem_num == 0){

		local::head = nullptr;	// of not element left, set head to be null. 

		LB::end = nullptr;
	}
	

}

///// @brief
///// Erase the elements in the sending list from the Hash table
//void Erase_elem_old(std::vector<int>::iterator& it)
//
//
//
//}
