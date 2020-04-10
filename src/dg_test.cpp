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
#include <mpi.h>

void Write_faces();


void Test(){

	if(mpi::rank == 2){

		for(auto& a : LB::proc_mapping_table){
		
			std::cout<< "rank "<<a.irank <<" gnum " << a.gnum << "\n" ;

		}

	}
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

void Write_faces_all(){

	for(int i = 0; i < mpi::num_proc; ++i){

		if(mpi::rank == i){
			Write_faces();

		}
		else{

			MPI_Barrier(MPI_COMM_WORLD);
		}

	}

}


void Write_faces(){
	

	Unit* temp =local::head;

		// generate the file name
	std::stringstream ss;
	ss << "../tests/facen" << std::setfill('0') << std::setw(5) << file_num1 << ".dat";
	std::string filename = 	ss.str();
	std::ofstream myfile; 	// stream class to write on files	

	if(mpi::rank == 0){
		myfile.open(filename, std::ios::trunc); // truncate the old file
	}
	else{
		myfile.open(filename, std::ios::app);

	}

	myfile<< "===============" << mpi::rank << "==============================="<< "\n";

	for(int k = 0; k < local::local_elem_num; ++k){
		int key = Get_key_fun(temp -> index[0], temp -> index[1], temp -> index[2]);
					
		myfile<< "---------------------------------------" << "\n";
	
		// coord
		myfile<< "elem: " << local::Hash_elem[key] -> index[0] <<
							 local::Hash_elem[key] -> index[1] <<
							 local::Hash_elem[key] -> index[2] << "\n";
		myfile<< "status "<< local::Hash_elem[key] -> status<< " c_position "<< 
					local::Hash_elem[key] -> child_position << "\n";
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
				
				// ref
				myfile << "ref_x: " << v.ref_x[0] << " "<< v.ref_x[1] << "\n";
				myfile << "ref_y: " << v.ref_y[0] << " "<< v.ref_y[1] << "\n";
			}
	
	
		}

		temp = temp -> next;
	}

		myfile.close();
		file_num1++;


	
}
