#include "dg_gen_status.h"
#include <cassert>
#include <cmath>
#include "dg_status_table.h"
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"
#include <iostream>	// test

// forward declaration----------------------------------------------
void Gen_status(int n, char* status);
//-------------------------------------------------------------------


/// @brief 
/// Generate the Hilbert status for all elements on root processor.
void Gen_status_all(){
	
	// allocate memory
	SortMesh::status = new char[SortMesh::num_of_element]{};

	Gen_status(grid::exp_x, SortMesh::status);
	


}


/// @brief 
/// Generate the Hilbert status. (Generate the status globally, on one processor.)
/// Restriction: square domain
/// @param n pow(2, n) : element number on each boundary in exponential form.
/// @param status element status array.
void Gen_status(int n, char* status){
	
	// make sure element number is valid.
	assert(n >= 0 && "n must >= 0");
	assert(mpi::rank == 0);
	

	// total element number
	const int total_elem = pow(2, 2 * n);
	
	const char four[]{'A', 'H', 'H', 'B', '\0'};
	
	// 3 situation, depending on the element number
	if(n > 1){	// more than 4 elem
		
		int* find = new int[n-1]{};	// record path

		for(int k = 0; k < total_elem; ++k){
			int div = total_elem / 4;
			unsigned int mid{};
			int s{};
			for(s = 0; s < n-1; ++s){
				mid = k % div;
				find[s] = k / div;
				div /= 4;	
			}
			--s;
			status[k] = four[find[s]];
			--s;

			for(s; s >= 0; --s){
				
				status[k] = Status_table(status[k], find[s]);

			}
			
			status[k] = Status_table(status[k], mid);

		}
	
		// free find[]
		delete[] find;		
		return;
	}
	else if(n == 1){	// 4 elem
		for(int k = 0; k < total_elem; ++k){
			status[k] = four[k];

		}
		return;
	}
	else{	// 1 elem

		status[0] = 'H';
		return;
	}



}
