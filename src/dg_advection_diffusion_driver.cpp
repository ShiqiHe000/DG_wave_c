#include <mpi.h>
#include "dg_param.h"
#include "dg_basis_storage.h"
#include "dg_create_hash.h"
#include "dg_construct_mpi_boundary.h"
#include "dg_init.h"
#include "dg_io.h"
#include "dg_adapt.h"
#include "dg_local_storage.h"
#include "dg_load_balancing.h"
#include "dg_derived_datatype.h"
#include "dg_step_by_RK3.h"
#include "dg_test.h"	// test
#include <iostream>	// test

/// @brief
/// Driver for DG approxiation. Algorithm 51. 
/// First, get DG basis parameters, such as collocation points and weights.
/// Then marches by each time step. Using explicit 3rd order Runge-Kutta methods.
void Driver_for_DG_approximation(){
	
	// construct basis
	Construct_basis_storage();

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
	
	Serial_io(tn);		

	Construct_data_type();

	// time integration
	for(int k = 0; k < dg_time::nt; ++k){

//		if(k == 0){
//			Adapt(k);
//			DG_init();	
//		}
		DG_step_by_RK3(tn, delta_t);

     //		Serial_io(tn);		
		Adapt(k);
//		Write_faces_all();
//     		Serial_io(tn);		

		// load_balancing----------------------------------------------	
//		Load_balancing(k);
		//-------------------------------------------------------------
//		Write_faces_all();

//		Simple_test(k);

      	//	Serial_io(tn);		
		tn = (k + 1) * delta_t;

	}

	Free_type();	// free up the derived datatype


}
