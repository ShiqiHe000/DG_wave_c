#include "dg_local_storage.h"
#include "dg_param.h"
#include "dg_unit.h"
#include "dg_user_defined.h"
#include <cmath>	// std::abs
#include "dg_verification.h"
#include <mpi.h>
#include <vector>
#include <iostream>	// test

/// @brief
/// Verfies your results. Now assume uniform grids. 
void Get_error(){

	result::L2_norm = std::vector<double>(dg_fun::num_of_equation);

	// temperary pointer
	Unit* temp = local::head;

	// traverse the linked list
	for(int k = 0; k < local::local_elem_num; ++k){
		
		// allocate
		int size_plane =(temp -> n + 1) * (temp -> m + 1);
		for(int i = 0; i < dg_fun::num_of_equation; ++i){
	
			result::error[i] = std::vector<double>(size_plane);
			result::exact[i] = std::vector<double>(size_plane);
		}

		// element size
		double del_x = (temp -> xcoords[1]) - (temp -> xcoords[0]); 
		double del_y = (temp -> ycoords[1]) - (temp -> ycoords[0]);  

		// wave ---------------------------------------------------------------------------------------		
//		Exact_solution_Gaussian(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
//					del_x, del_y, result::exact, dg_time::t_total);
		// --------------------------------------------------------------------------------------------

		// test ---------------------------------------------------------------------------------------
		Exact_solution_sin(temp -> n, temp -> m, temp -> xcoords[0], temp -> ycoords[0], 
					del_x, del_y, result::exact, dg_time::t_total);
		//----------------------------------------------------------------------------------------------

		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
	
			for(int i = 0; i < size_plane; ++i){
		
				result::error[equ][i] = std::abs(result::exact[equ][i] - temp -> solution[equ][i]);

				result::L2_norm[equ] += result::error[equ][i] * result::error[equ][i];
				
			}

		}

		result::error.clear();
		result::exact.clear();

		temp = temp -> next;	

	}

	// sum up L2_norm on processor 0
	std::vector<double> L2_recv(dg_fun::num_of_equation);
	MPI_Reduce(&result::L2_norm[0], &L2_recv[0], dg_fun::num_of_equation, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(mpi::rank == 0){
		for(int equ = 0; equ < dg_fun::num_of_equation; ++equ){
			result::L2_norm[equ] = sqrt(L2_recv[equ]);
			std::cout<< result::L2_norm[equ] << "\n";
		}
	}

}

