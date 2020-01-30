#include "dg_test.h"
#include <iostream>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_param.h"
#include "dg_cantor_pairing.h"

void Write_faces(int key);


void Test(){

//	Unit* temp = local::head;
//	if(mpi::rank == 0){
//		for(int k = 0; k < local::local_elem_num; ++k){
//			
//	
//			for(int i=0; i < 4; ++i){
//	
//				std::cout<<"elem" << k << "\n";
//				std::cout<< i << " " <<temp -> facen[i][0].face_type << "\n";
//				std::cout<< i << " " <<temp -> facen[i][0].hlevel << "\n";
//				std::cout<< i << " " <<temp -> facen[i][0].porder << "\n";
//				std::cout<< i << " " <<temp -> facen[i][0].key << "\n";
//				std::cout<< i << " " << temp -> facen[i][0].rank << "\n";
//			}
//			temp = temp -> next;
//		}
//	}

//	for(int i = 0; i <= mpi::num_proc; ++i){
//
//		std::cout<< mpi::rank << " " << local::elem_range[i] << "\n";
//
//	}

	Unit* temp =local::head;
	if(mpi::rank == 2){

		for(int k = 0; k < local::local_elem_num; ++k){
			int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
						

			Write_faces(key);
			
			temp = temp -> next;
		}

	}

}

void Write_faces(int key){
	
	std::cout<< "---------------------------------------" << "\n";

	// coord
	std::cout<< "elem: " << local::Hash_elem[key] -> index[0] <<
						 local::Hash_elem[key] -> index[1] <<
						 local::Hash_elem[key] -> index[2] << "\n";
	// write 4 faces
	for(int i = 0; i < 4; ++i){

		if(i == 0){
			std::cout<< "South===================" << "\n";
		}
		else if(i == 1){
			std::cout<< "North===================" << "\n";
		}
		else if(i == 2){

			std::cout<< "West===================="<< "\n";
		}
		else{
			std::cout<< "East===================="<< "\n";

		}

		
		for(auto& v : local::Hash_elem[key] -> facen[i]){
				

			//face type
			std::cout<< "face_type: " << v.face_type << "\n";

			// hlevel
			std::cout<< "hlevel: " << v.hlevel << "\n";

			// porder
			std::cout<< "porder: " << v.porderx << " " << v.pordery << "\n";

			// key
			std::cout << "key: "<< v.key << "\n";

			// rank
			std::cout << "rank: " << v.rank << "\n";

		}


	}

}
