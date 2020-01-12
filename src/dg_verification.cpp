#include "dg_local_storage.h"
#include "dg_param.h"
#include "dg_unit.h"
#include "dg_user_defined.h"
#include <cmath>	// std::abs
#include "dg_verification.h"
#include <mpi.h>
#include <iostream>	// test

/// @brief
/// Verfies your results. Now assume uniform grids. 
void Get_error(){

	// allocate
	int size_err =(grid::nmin + 1) * (grid::nmin + 1) * dg_fun::num_of_equation;
	int size_plane =(grid::nmin + 1) * (grid::nmin + 1);
	result::error = new double[size_err]{};
	result::exact = new double[size_err]{};
	result::L2_norm = new double[dg_fun::num_of_equation]{};

	// temperary pointer
	Unit* temp = local::head;

	// traverse the linked list
	for(int k = 0; k < local::local_elem_num; ++k){
		
		// element size
		double del_x = (temp -> xcoords[1]) - (temp -> xcoords[0]); 
		double del_y = (temp -> ycoords[1]) - (temp -> ycoords[0]);  
		
		Exact_solution_Gaussian(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
					del_x, del_y, result::exact, dg_time::t_total);


		int nodei{};
		int normi{};
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
			for(int i = 0; i < size_plane; ++i){
		
				result::error[nodei] = std::abs(result::exact[nodei] - temp -> solution[nodei]);

				result::L2_norm[normi] += result::error[nodei] * result::error[nodei];
				
				++nodei;
			}

			++normi;
		}

		temp = temp -> next;	

	}

	MPI_Reduce(result::L2_norm, L2_recv, dg_fun::num_of_equation, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(mpi::rank == 0){
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
			result::L2_norm[equ] = sqrt(L2_recv[equ]);
			std::cout<< L2_norm[equ] << "\n";
		}
	}


	delete[] result::error;
	delete[] result::exact;
	delete[] result::L2_norm;
}

