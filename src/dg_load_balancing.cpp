#include "dg_load_balancing.h"
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_param.h"
#include <cmath>
#include <mpi.h>
#include <cassert>
#include <iostream> // test

// forward declaration -----------------------------------------
double Elem_load(int porder);
// ------------------------------------------------------------

/// @brief Calculate the sum of the local computational load.
void Sum_up_local_load(){

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
	int proc_pre = -1;
	for(int k = 0; k < local::local_elem_num; ++k){
		
		LB::lprefix_load[k] += exscan_sum;

		LB::pmapping[k] = std::floor((LB::lprefix_load[k] - 0.01 * load_avg) / load_avg);

		if(LB::pmapping[k] != proc_pre){

			LB::proc_mapping_table.push_back({LB::pmapping[k], k + elem_accum});
			
			proc_pre = LB::pmapping[k];
		}	

		
//--------test--------------------------------------------------------------------------		
		assert(LB::pmapping[k] > 0 && "preocessor mapping is small than 0.");
//--------------------------------------------------------------------------------------

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

		int first_rank = LB::proc_mapping_table.begin(),irank;

		if(first_rank == pre_rank){	// if equal than erase the first column

			LB::proc_mapping_table.erase(LB::proc_mapping_table.begin());
		}

	}
	
	// wait
	if(mpi::rank != (mpi::num_proc - 1)){
	
		MPI_Status status2;
		MPI_Wait(&request, &status2);
	}


	

}

/// @brief 
/// Element computational load. The load on each element due to fluid computations is O(N**4),
/// where N is the number of grid points along one direction. The load is normalized between 0 and 1. 
/// @param porder input the polynomial order of current element. (Now we assume polynomial
/// order are identical on two direction.)
double Elem_load(int porder){

	static const double load_min = std::pow((double)(grid::nmin + 1), 4);
	static const double load_max = std::pow((double)(grid::nmax + 1), 4);

	double load = (std::pow(porder + 1, 4) - load_min) / (load_max - load_min);

	return load;

}
