#include "dg_gen_status.h"
#include <cassert>
#include <cmath>
#include "dg_status_table.h"
#include "dg_param.h"
#include "dg_nodal_2d_storage.h"
#include <vector>
#include <iostream>	// test

// forward declaration----------------------------------------------
void Gen_status(int n, std::vector<char>& status);
//-------------------------------------------------------------------


/// @brief 
/// Generate the Hilbert status for all elements on root processor.
void Gen_status_all(){
	
	// allocate memory
	SortMesh::status = std::vector<char>(SortMesh::num_of_element);

	Gen_status(grid::exp_x, SortMesh::status);

}


/// @brief 
/// Generate the Hilbert status. (Generate the status globally, on one processor.)
/// Restriction: square domain
/// @param n pow(2, n) : element number on each boundary in exponential form.
/// @param status element status array.
void Gen_status(int n, std::vector<char>& status){
	
	// make sure element number is valid.
	assert(n >= 0 && "n must >= 0");
	assert(mpi::rank == 0);
	

	// total element number
	std::vector<char> pre{'A', 'H', 'H', 'B'};

	if(n > 1){
		
		for(int k = 2; k <= n; ++k){

			int now_num = (int)std::pow(2, 2 * k);
			std::vector<char> inter(now_num);

			auto it = inter.begin();

			for(auto& v: pre){

				for(int i = 0; i < 4; ++i){
					*(it) = Status_table(v, i);
					++it;
				}
			}
			
			pre.clear();
			pre = inter;

		}
		status = pre;
		return;

	}
	else if(n == 1){

		status = pre;

		return;
	}
	else{	// n == 0

		status[0] = 'H';
		return;

	}



}
