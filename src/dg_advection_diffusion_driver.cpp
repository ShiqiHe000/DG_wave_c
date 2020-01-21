#include <mpi.h>
#include "dg_param.h"
#include "dg_constructor.h"
#include "dg_create_hash.h"
#include "dg_construct_mpi_boundary.h"
#include "dg_init.h"
#include "dg_io.h"
#include "dg_adapt.h"

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
	
	// h-refinement
	h_refinement();
	

	// time integration
//	for(int k = 0; k < dg_time::nt; ++k){
//	
//		Serial_io(tn);		
//
//		tn = (k + 1) * delta_t;
//
//	}


}
