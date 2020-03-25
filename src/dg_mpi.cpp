#include <mpi.h>
#include "dg_param.h"	// mpi variables
#include "dg_mpi.h"

/// Initialize mpi environment. 
/// @param argc pointer to number of argument vector
/// @param argv array of pointer to argument vector
void Start_mpi(int argc, char** argv){
	
	// Initialize
	MPI_Init(&argc, &argv);
	
	// MPI env
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi::num_proc); 


}
