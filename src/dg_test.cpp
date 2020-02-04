#include "dg_test.h"
#include <iostream>
#include "dg_unit.h"
#include "dg_local_storage.h"
#include "dg_param.h"
#include "dg_cantor_pairing.h"
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <iomanip>

void Write_faces( int nt);


void Test(int nt){

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


//	if(mpi::rank == 0){
//
//		std::unordered_map<int, Unit*>::const_iterator got = local::Hash_elem.find(4);
//
//
//		if(got == local::Hash_elem.end()){
//
//			std::cout<< "not find" << "\n";
//
//		}
//
//	}

//	Unit* temp =local::head;
//	if(mpi::rank == 0){
//
//		for(int k = 0; k < local::local_elem_num; ++k){
//			int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
//						
//
//			Write_faces(key, nt);
//			
//			temp = temp -> next;
//		}
//
//	}

}

int file_num1 = 1;

void Write_faces(int nt){
	

	Unit* temp =local::head;
	if(mpi::rank == 0){

		// generate the file name
		std::stringstream ss;
		ss << "../tests/facen" << std::setfill('0') << std::setw(5) << file_num1 << ".dat";
		std::string filename = 	ss.str();
		std::ofstream myfile; 	// stream class to write on files	

		myfile.open(filename, std::ios::trunc); // truncate the old file


		for(int k = 0; k < local::local_elem_num; ++k){
			int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
						
			myfile<< "---------------------------------------" << "\n";
		
			// coord
			myfile<< "elem: " << local::Hash_elem[key] -> index[0] <<
								 local::Hash_elem[key] -> index[1] <<
								 local::Hash_elem[key] -> index[2] << "\n";
			// write 4 faces
			for(int i = 0; i < 4; ++i){
		
				if(i == 0){
					myfile<< "South===================" << "\n";
				}
				else if(i == 1){
					myfile<< "North===================" << "\n";
				}
				else if(i == 2){
		
					myfile<< "West===================="<< "\n";
				}
				else{
					myfile<< "East===================="<< "\n";
		
				}
		
				
				for(auto& v : local::Hash_elem[key] -> facen[i]){
						
		
					//face type
					myfile<< "face_type: " << v.face_type << "\n";
		
					// hlevel
					myfile<< "hlevel: " << v.hlevel << "\n";
		
					// porder
					myfile<< "porder: " << v.porderx << " " << v.pordery << "\n";
		
					// key
					myfile << "key: "<< v.key << "\n";
					
					// rank
					myfile << "rank: " << v.rank << "\n";
		
				}
		
		
			}

			temp = temp -> next;
		}

		myfile.close();
		file_num1++;


	}
}
