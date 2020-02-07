#include "dg_load_balancing.h"
#include "dg_unit.h"
#include "dg_local_storage.h"


/// @brief Calculate the sum of the local computational load.
void Sum_up_local_load(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){

		
		
		temp = temp -> next;
	}	

}
