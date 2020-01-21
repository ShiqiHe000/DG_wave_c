#include "dg_test.h"
#include <iostream>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_param.h"

void Test(){

	Unit* temp = local::head;
	if(mpi::rank == 1){
		for(int k = 0; k < local::local_elem_num; ++k){
			
	
			for(int i=0; i < 4; ++i){
	
				std::cout<<"elem" << k << "\n";
				std::cout<< i << " " <<temp -> facen[i][0].face_type << "\n";
				std::cout<< i << " " <<temp -> facen[i][0].hlevel << "\n";
				std::cout<< i << " " <<temp -> facen[i][0].porder << "\n";
				std::cout<< i << " " <<temp -> facen[i][0].key << "\n";
	
			}
			temp = temp -> next;
		}
	}
}
