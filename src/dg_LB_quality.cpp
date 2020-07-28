#include <mpi.h>
#include "dg_param.h"
#include "dg_LB_quality.h"
#include "dg_element_load.h"
#include <vector>
#include "dg_local_storage.h"
#include <fstream>	// io
#include <iomanip>	// std::setw
#include <sstream>
#include <iostream>	// test


/// @brief
/// Evaluate the workload balancing efficiecy of the current time. 
/// @param t current time step 
void LB_efficiency(double t){

	Unit* temp = local::head;

	std::vector<double> lprefix_load(local::local_elem_num);

	// local prefix sum of load
	lprefix_load[0] = Elem_load(temp -> n);
	temp = temp -> next;
	for(int k = 1; k < local::local_elem_num; ++k){
		lprefix_load[k] = Elem_load(temp -> n) + lprefix_load[k - 1];
		
		temp = temp -> next;

	}	
	
	double local_load_sum = lprefix_load.back();	// local computational load sum
	double exscan_sum{};	// the load of former processor
//if(mpi::rank == 0){
//
//	std::cout<< "LB quality: local load sum "<< local_load_sum << "\n";
//}
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
	

	// rank 0 get the max local workload 
	double max_local_load{};
	MPI_Reduce(&local_load_sum, &max_local_load, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
//	MPI_Allreduce(&local_load_sum, &max_local_load, 1,
//                  	MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

//if(local_load_sum == max_local_load){
//
//	std::cout<< "time "<< t << " rank " << mpi::rank << "\n";
//}


	// compute the efficiency
	if(mpi::rank == 0){
	
		double eff = load_avg / max_local_load * 100.0;

		double opt_balance{};

		if(LB::opt_bottleneck > 0){
			opt_balance = load_avg / LB::opt_bottleneck * 100.0;

		}

		std::stringstream ss;
		ss << fileinfo::eff_filename << ".csv";
		std::string filename = 	ss.str();
		std::ofstream eff_file;


		if(t == 0){	// first time step
			eff_file.open(filename, std::ios::trunc);

			eff_file << "Time(s) Efficiency(%) Optimal_balance\n";
		}
		else{

			eff_file.open(filename, std::ios::app);

		}

		eff_file << t << " " << eff << " " << opt_balance << "\n";

	
		eff_file.close();
			
	}

	LB::opt_bottleneck = 0.0;
}
