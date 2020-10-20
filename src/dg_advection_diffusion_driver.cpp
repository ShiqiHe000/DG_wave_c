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
#include "dg_CFL.h"
#include "dg_LB_quality.h"	// LB quality
#include "dg_test.h"	// test
#include "dg_reinit.h"	// reinit the solution
#include <iostream>	// test

/// @brief
/// Driver for DG approxiation. Algorithm 51. 
/// First, get DG basis parameters, such as collocation points and weights.
/// Then marches by each time step. Using explicit 3rd order Runge-Kutta methods.
void Driver_for_DG_approximation(){
	
	// check CFL condition
	Check_CFL();

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
	DG_init();	// init with exact solution
//	DG_init2();	// init with exact solution
//	DG_init_new();	// init with other value;
	
	Serial_io(tn);		

//	LB_efficiency_write(tn);

	Construct_data_type();
	
	// time integration
	for(int k = 0; k < dg_time::nt; ++k){

		DG_step_by_RK3(tn, delta_t);
		
		tn = (k + 1) * delta_t;

//		if(k == 0){
//			DG_reinit(tn);
//		}
	

		// output control
		if(dg_io::io){
			if((k + 1) % dg_io::output_frequency == 0){
	     			Serial_io(tn);		
			}
		}

		if(dg_refine::adapt){	// hp-refinement

			if((k + 1) % dg_refine::refine_frequency == 0){
				// hp-adaptive --------------------------------------------
				Adapt(k);
//				LB_efficiency_write(tn);
				// --------------------------------------------------------
			
				// reinit solutions---------------------------------------
				if(k == 0){
					DG_reinit(tn);
				}
				// --------------------------------------------------------


//     				Serial_io(tn);		
				if(dg_refine::load_balancing){	// repartitioning

					
					// load_balancing----------------------------------------------	
					Load_balancing(k);
					
					LB_efficiency_evaluate();
				//	Write_faces_all();
					//-------------------------------------------------------------

					int LB_times{};

					while(!LB::high_eff){

						++LB_times;
				
						Load_balancing(k);
						
						LB_efficiency_evaluate();

						if(LB_times > 3){	// avoid deadlock

							break;	
						}

					}

					LB_set_back();

//					LB_efficiency_write(tn);
//     					Serial_io(tn);		
				}
//				Write_faces_all();
			}
		}

	}

	Free_type();	// free up the derived datatype


}
