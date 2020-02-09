#include "dg_load_balancing.h"
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_param.h"
#include <cmath>
#include <mpi.h>
#include "dg_cantor_pairing.h"
#include <iostream> // test

// forward declaration -----------------------------------------
double Elem_load(int porder);
// ------------------------------------------------------------

// global variable----------------------------------------------
struct sending_envelope Send;	// record what to send
//--------------------------------------------------------------

/// @brief Calculate the sum of the local computational load.
void Build_mapping_table(){

	Unit* temp = local::head;

	LB::lprefix_load = std::vector<double> (local::local_elem_num);
	LB::pmapping = std::vector<int> (local::local_elem_num);

	// local prefix sum of load
	LB::lprefix_load[0] = Elem_load(temp -> n);
	temp = temp -> next;
	for(int k = 1; k < local::local_elem_num; ++k){
		LB::lprefix_load[k] = Elem_load(temp -> n) + LB::lprefix_load[k - 1];
		
		temp = temp -> next;

	}	
	
	double local_load_sum = LB::lprefix_load.back();	// local computational load sum
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
	int elem_accum{};
	MPI_Exscan(&local::local_elem_num, &elem_accum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	// form processor mapping table
	temp = local::head;
	int proc_pre = - 1;
	for(int k = 0; k < local::local_elem_num; ++k){
		
		LB::lprefix_load[k] += exscan_sum;

		LB::pmapping[k] = std::floor((LB::lprefix_load[k] - 0.01 * load_avg) / load_avg);
		
		// form partial mapping table
		if(LB::pmapping[k] != proc_pre){

			LB::proc_mapping_table.push_back({LB::pmapping[k], k + elem_accum});
			
			proc_pre = LB::pmapping[k];
		}	
		
		// form sending_envelope
		if(LB::pmapping[k] < mpi::rank){	// need to be moved to the former proc

			int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

			Send.pre.push_back(key);
		}
		else if(LB::pmapping[k] > mpi::rank){

			int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);

			Send.next.push_back(key);
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

	// rebild the mapping table
	LB::proc_mapping_table.clear();	// clear all the elements
	for(int i = 0; i < mpi::num_proc; ++i){
		
		LB::proc_mapping_table.push_back({recv_buff[2 * i], recv_buff[2 * i + 1]});

	}


}

/// @brief
/// After built the complete mapping table, now we decide how to reallocate the elements
/// @param elem_accum accumulated element of former processors.
void Reallocate_elem(int elem_accum){

	int first = elem_accum; 	// first elem global number
	int last = first + local::local_elem_number - 1;	// last elem global number

	if(mpi::rank == 0){	// first proc
		
		if(LB::proc_mapping_table[1].gnum <= last){	// should send

			int num_send = last - LB::proc_mapping_table[1].gnum + 1;	// num of element to be sent

		}
		else if((LB::proc_mapping_table[1].gnum - 1) > last){	// recv from r1


		}

		// otherwiase no send or recv

	}
	else if(mpi::rank == mpi::num_proc - 1){	// last proc
		
		if(LB::proc_mapping_table.back().gnum < start){	// should recv


		}
		else if(LB::proc_mapping_table.back().gnum > start){	// send

			int num_send = LB::proc_mapping_table.back().gnum - start;

		}

	}
	else{	// proc in between 

		// with former proc
		if(LB::proc_mapping_table[mpi::rank] < start){	// recv


		}
		else if(LB::proc_mapping_table[mpi::rank] > start){	// send

			int num_send = LB::proc_mapping_table[mpi::rank] - start;

		}

		// with latter proc
		if(LB::proc_mapping_table[mpi::rank + 1] <= last){	// send   

			int num_send = LB::proc_mapping_table[mpi::rank + 1] - last + 1;

		}
		else if(LB::proc_mapping_table[mpi::rank + 1] > (last - 1)){	// recv


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
