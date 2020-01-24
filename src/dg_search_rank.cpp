#include "dg_search_rank.h"
#include "dg_local_storage.h"
#include "dg_param.h"
#include "dg_hilbert_curve.h"

/// @brief
/// Search neighbour element's rank based on it's integer
/// coordinates. Only works for uniform grid. 
/// @param i neighbour's x coordinates
/// @param j neighbour's y coordinates
int Target_rank(int i, int j){

	// neighbour's Hilbert number
	int d = xy2d(grid::exp_x, j, i);

	int target_rank = mpi::rank;

	if(d > local::elem_range[mpi::rank + 1]){
		++target_rank;

		for(int k = mpi::rank + 2; k <= mpi::num_proc; ++k){
			if(d > local::elem_range[k]){
				++target_rank;
			}
			else{
				return target_rank;
			}

		}

	}
	else{

		--target_rank;

		for(int k = mpi::rank - 1; k >= 1; --k){

			if(d <= local::elem_range[k]){

				--target_rank;

			}
			else{
				return target_rank;
			}

		}

	}

}

