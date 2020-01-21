//#include "dg_adapt.h"
#include "dg_unit.h"
#include "dg_local_storage.h"
#include <cstdlib>	// random number

/// @brief
/// Random h-refinement scheme. Each element has 30% chance to split.
void h_refinement(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		// generate random number
		int rand_num = rand() % 10 + 1;	// random number between [1, 10]

		if(rand_num <= 3){	// refine

			
		}


		temp = temp -> next;
	}

}
