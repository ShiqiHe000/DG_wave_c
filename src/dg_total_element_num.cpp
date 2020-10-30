#include <iostream>
#include <mpi.h>
#include "dg_local_storage.h"
#include "dg_total_element_num.h"
#include "dg_param.h"

/// @brief
/// Root node 0 get the total number of element. 
/// @param t Current time. 
void Total_element_num(double t){

	int total{};

	MPI_Reduce(&local::local_elem_num, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if(mpi::rank == 0){

		std::cout<< "Total num of elements at time " << t << " is " << total << "\n";
	}
}
