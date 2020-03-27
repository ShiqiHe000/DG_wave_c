#include <vector>
#include "dg_local_storage.h"


void Numerical_flux_x(){

	Unit* temp = local::head;

	for(int k = 0; k < local::local_elem_num; ++k){
		
		// compute numerical flux on the south interface
		for(auto it_face = it -> facen[0].begin(); it_face != it -> facen[0].end(); ++it_face){

			if(it_face -> face_typy == 'L'){	// local neighbour

			}
			else if(it_face -> face_type == 'B'){	// phsical boundary

				int size = dg_fun::num_of_equation * (temp -> m + 1);
				std::vector<double> solution_ext(size);
			}
			else{	// mpi boundary


			}
		}


		temp = temp -> next;
	}
}
