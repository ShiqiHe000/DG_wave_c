#include <mpi.h>
#include "dg_param.h"
#include "dg_constructor.h"
#include "dg_create_hash.h"
#include "dg_construct_mpi_boundary.h"
#include "dg_init.h"

/// @brief
/// Driver for DG approxiation. Algorithm 51. 
/// First, get DG basis parameters, such as collocation points and weights.
/// Then marches by each time step. Using explicit 3rd order Runge-Kutta methods.
void Driver_for_DG_approximation(){

	// construct basis
	Construct_basis();

	// create hash table
	Create_hash();

	// construct mpi boundaries and physical boundaries
	MPI_boundary_construct();

	// time step
	double delta_t = dg_time::t_total / dg_time::nt;

	// current time
	double tn{};

	// Initialization
	DG_init();	
	



}
