#include <mpi.h>
#include "dg_param.h"
#include "dg_LB_quality.h"
#include "dg_element_load.h"


/// @brief
/// Evaluate the workload balancing efficiecy of the current time. 
/// @param t current time step 
void LB_efficiency(double t){

//	Unit* temp = local::head;
//
//	std::vector<double> lprefix_load(local::local_elem_num);
//
//	// local prefix sum of load
//	lprefix_load[0] = Elem_load(temp -> n);
//	temp = temp -> next;
//	for(int k = 1; k < local::local_elem_num; ++k){
//		lprefix_load[k] = Elem_load(temp -> n) + lprefix_load[k - 1];
//		
//		temp = temp -> next;
//
//	}	
//	
//	double local_load_sum = lprefix_load.back();	// local computational load sum
//	double exscan_sum{};	// the load of former processor
//	
//	// Global prefix sum of load
//	MPI_Exscan(&local_load_sum, &exscan_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//	
//	// calculate for the average load
//	double load_avg{};
//	if(mpi::rank == (mpi::num_proc - 1)){ // last proc does the job
//
//		double load_tol = exscan_sum + local_load_sum;
//
//		load_avg = load_tol / mpi::num_proc;
//
//	}
//	
//	// broadcast average load
//	MPI_Bcast(&load_avg, 1, MPI_DOUBLE, mpi::num_proc - 1, MPI_COMM_WORLD);
//	
//	// Global element number
//	MPI_Exscan(&local::local_elem_num, &LB::elem_accum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//

}
